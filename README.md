# WMM2025C

World Magnetic Model (WMM-2025) C Implementation

## Overview

This project provides a C library and test executable for calculating geomagnetic field elements using the World Magnetic Model (WMM-2025). It supports reading WMM COF files, computing magnetic elements (declination, inclination, total intensity, grid variation), Cartesian components, annual change (secular variation), and representative uncertainties.

- **Reference:** Based on NOAA/NCEI WMM2020 Legacy C code

## Features

- Reads WMM COF files (model coefficients)
- Computes geomagnetic elements: Declination (D), Inclination (I), Total Intensity (F), Grid Variation (GV)
- Computes Cartesian components: X (North), Y (East), Z (Down), H (Horizontal)
- Calculates annual change (secular variation) for all elements
- Includes a test executable for validation

## Build Instructions

### Prerequisites

- CMake >= 3.15
- Tested compiler: MSVC 2022

### Build Steps

1. Clone the repository and enter the project directory:

    ```sh
    git clone <repo-url>
    cd wmm2025C
    ```

2. Configure and build (default: Debug):

    ```sh
    cmake -S . -B build
    cmake --build build
    ```

    To build in Release mode:

    ```sh
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build
    ```

3. The test executable (`test_wmm_mag`) will be built in the `build` directory.

### Library Options

- By default, builds a **static library**.
- To build a **shared library**:

    ```sh
    cmake -S . -B build -DWMM_MAG_BUILD_SHARED=ON
    ```

## Usage

- **Library:** Link your C project against `wmm_mag` and include `wmm_mag.h`.
- **Test Executable:** Run `test_wmm_mag` to validate model calculations.

## Data Files

- Place WMM COF files (e.g., `WMM2025.COF`) and test vectors (`WMM2025_TEST_VALUES.txt`) in the project root. These files are automatically copied to the build/run directory.

## File Structure

- `wmm_mag.c` / `wmm_mag.h` : Main library source and API
- `wmm_mag_struct.h`        : Data structures for model/state
- `test_wmm_mag.c`          : Test executable source
- `CMakeLists.txt`          : Build configuration
- `*.COF`                   : WMM model coefficient files
- `WMM2025_TEST_VALUES.txt` : Test vectors


## Example

```C
const char *cof = "WMM2025.COF"; /* NOAA/NCEI COF file */

// WMM initialization (reset state, load COF, internal init)
wmm_init(cof);

// Get WMM Version
double epoch;
char name[32];
char date[16];
wmm_version(&epoch, name, sizeof(name), date, sizeof(date));

// Calculate WMM field at a specific location and time
double decl, incl, ti1, gv;
double x1, y1, z1, h1;

wmm_point(alt, dlat, dlon, time,
          &decl, &incl, &ti1, &gv,
          &x1, &y1, &z1, &h1);
```

## References

- [WMM2025 Technical Report](https://www.ngdc.noaa.gov/geomag/WMM/)
- [NOAA/NCEI WMM2020 Legacy C code](https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2020/WMM2020LegacyC.zip)

    ```
    NOTICE — Use of U.S. Government Work (WMM)

    This product includes or is derived from the World Magnetic Model (WMM)
    information and source code produced by the U.S. Government (NOAA/NCEI).
    The WMM source code is in the public domain and not licensed or under copyright.
    The information and software may be used freely by the public.

    Pursuant to 17 U.S.C. §403, third parties producing copyrighted works that
    consist predominantly of U.S. Government material must provide notice
    identifying the U.S. Government material incorporated and stating that such
    material is not subject to copyright protection.
    ```