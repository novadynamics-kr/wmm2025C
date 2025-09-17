cmake_minimum_required(VERSION 3.15)
project(WMM2025C)

# =========================================================
# Default build type (change with -DCMAKE_BUILD_TYPE=Release if desired)
# =========================================================
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
endif()
message(STATUS "Active build configuration: ${CMAKE_BUILD_TYPE}")

# Common compile options
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS_DEBUG   "-Wall -Wextra -W -g -O0")
set(CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -W -O2")

if(MSVC)
  set(CMAKE_C_FLAGS_DEBUG   "/W4 /Zi /Od /wd4305")
  set(CMAKE_C_FLAGS_RELEASE "/W4 /O2 /wd4305")
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
  add_compile_options(/utf-8)
  message(STATUS "Compiler : MSVC")
endif()

# =========================================================
# Option: Choose library type (default: static)
#   Build shared (dynamic) library with -DWMM_MAG_BUILD_SHARED=ON
# =========================================================
option(WMM_MAG_BUILD_SHARED "Build wmm_mag as a shared library" ON)

# =========================================================
# Library target: wmm_mag (static or shared)
# =========================================================
set(WMM_MAG_SOURCES
    wmm_mag.c
)
set(WMM_MAG_PUBLIC_HEADERS
    wmm_mag.h
    wmm_mag_struct.h
)

if(WMM_MAG_BUILD_SHARED)
  add_library(wmm_mag SHARED ${WMM_MAG_SOURCES})
  # DLL 빌드 시 WMM_MAG_EXPORTS 정의
  target_compile_definitions(wmm_mag PRIVATE WMM_MAG_EXPORTS)
else()
  add_library(wmm_mag STATIC ${WMM_MAG_SOURCES})
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(wmm_mag PRIVATE m)
  message(STATUS "Compiler : GCC")
endif()

# Set public interface include path for other projects
target_include_directories(wmm_mag PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Enable PIC for static library (useful for some platforms/linkers)
set_target_properties(wmm_mag PROPERTIES POSITION_INDEPENDENT_CODE ON)

# =========================================================
# Test executable: test_wmm_mag
#   - Links to wmm_mag library
# =========================================================
add_executable(test_wmm_mag test_wmm_mag.c)
target_link_libraries(test_wmm_mag PRIVATE wmm_mag)
target_include_directories(test_wmm_mag PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# =========================================================
# Runtime data (.COF, test vectors) copy
#   - Copy COF/text files to build/run folder
# =========================================================
file(GLOB COF_FILES "${CMAKE_SOURCE_DIR}/*.COF")
file(GLOB TXT_FILES "${CMAKE_SOURCE_DIR}/WMM2025_TEST_VALUES.txt")

set(MY_RUNTIME_TARGETS test_wmm_mag)  # List only executables

foreach(tgt ${MY_RUNTIME_TARGETS})
  foreach(file ${COF_FILES})
    add_custom_command(TARGET ${tgt} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${file} $<TARGET_FILE_DIR:${tgt}>
    )
  endforeach()
  foreach(file ${TXT_FILES})
    add_custom_command(TARGET ${tgt} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${file} $<TARGET_FILE_DIR:${tgt}>
    )
  endforeach()
endforeach()
