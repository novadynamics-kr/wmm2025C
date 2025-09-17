/*======================================================================*
 *  File        : wmm_mag.h
 *  Brief       : World Magnetic Model (WMM-2025) public API (C)
 *  Author      : Nova Dynamics (novadynamics.kr@gmail.com)
 *  License     : MIT
 *
 *  Notes
 *  - This header exposes a minimal, stable API to evaluate WMM fields.
 *  - Units:
 *      * latitude/longitude: degrees (geodetic; East/North positive)
 *      * altitude          : kilometers above WGS84 ellipsoid
 *      * time              : decimal year (e.g., 2025.0)
 *      * field outputs     : nT; angles in degrees
 *======================================================================*/

#ifndef WMM_MAG_H
#define WMM_MAG_H

#include <stddef.h> // for size_t

#ifdef _WIN32
  #ifdef WMM_MAG_EXPORTS
    #define WMM_MAG_API __declspec(dllexport)
  #else
    #define WMM_MAG_API __declspec(dllimport)
  #endif
#else
  #define WMM_MAG_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*======================================================================*
 *  Initialization / Model loading
 *======================================================================*/

/**
 * @brief Initialize the global WMM state from a COF file.
 *
 * @param[in]     cof     Path to WMM COF file (e.g., "WMM2025.COF").
 *
 * @details
 * This is a lightweight initializer that resets internal state,
 * loads Gauss coefficients (and secular-variation terms) from the COF,
 * and prepares Schmidt-semi-normalized tables used during synthesis.
 */
WMM_MAG_API void wmm_init(const char *cof);

/**
 * @brief Return the current WMM model metadata (epoch, model name, release date).
 *
 * Copies the model information that was loaded from the COF file during
 * initialization into caller-provided outputs. All output arguments are
 * optional—pass `NULL` (and/or a capacity of 0 for strings) to skip any field.
 *
 * @pre Call ::wmm_init() successfully before using this function, so that
 *      the global metadata (epoch, model name, release date) is populated.
 *
 * @param[out] epoch
 *     Pointer to a double that receives the model epoch as a decimal year
 *     (e.g., `2025.0`). May be `NULL`.
 * @param[out] model_name
 *     Destination buffer for the model name string (UTF-8), e.g., `"WMM-2025"`.
 *     May be `NULL` if you do not need the name.
 * @param[in]  model_name_cap
 *     Capacity of @p model_name in bytes (including space for the terminator).
 *     Ignored if @p model_name is `NULL`. If > 0, the result is always
 *     NUL-terminated.
 * @param[out] release_date
 *     Destination buffer for the release date string in `"MM/DD/YYYY"` format.
 *     May be `NULL` if you do not need the date.
 * @param[in]  release_date_cap
 *     Capacity of @p release_date in bytes (including space for the terminator).
 *     Ignored if @p release_date is `NULL`. If > 0, the result is always
 *     NUL-terminated.
 *
 * @note Typical buffer sizes used by this library are 32 bytes for the model
 *       name and 16 bytes for the release date. Larger buffers are fine.
 * @warning This function reads from a global metadata object initialized by
 *          ::wmm_init(); concurrent reinitialization from another thread is not
 *          safe unless externally synchronized.
 *
 * @par Example
 * @code
 *   double epoch;
 *   char name[32];
 *   char date[16];
 *   wmm_init("WMM2025.COF");
 *   wmm_version(&epoch, name, sizeof(name), date, sizeof(date));
 *   // epoch -> 2025.0, name -> "WMM-2025", date -> e.g., "12/15/2024"
 * @endcode
 */
WMM_MAG_API void wmm_version(double *epoch,
                 char *model_name, size_t model_name_cap,
                 char *release_date, size_t release_date_cap);


/*======================================================================*
 *  Time helper
 *======================================================================*/

/**
 * @brief Convert calendar date/time (UTC) to decimal year.
 *
 * @param[in]  year    e.g., 2025
 * @param[in]  month   1..12 (if 0, returns @p year as-is; legacy mode)
 * @param[in]  day     1..MonthDays
 * @param[in]  hour    0..23
 * @param[in]  minute  0..59
 * @param[out] decimal_year  decimal year (e.g., 2025.500…)
 * @return 1 on success, 0 on invalid inputs.
 *
 * @note Leap year handled per Gregorian rules. Seconds are ignored.
 *
 */
WMM_MAG_API int dec_year(int year, int month, int day, int hour, int minute, double *decimal_year);

/*======================================================================*
 *  Core evaluation
 *======================================================================*/

/**
 * @brief Compute main geomagnetic elements and NED components at a point.
 *
 * @param[in]  altitude_km     Altitude above WGS84 ellipsoid (km).
 * @param[in]  latitude_deg    Geodetic latitude (deg).
 * @param[in]  longitude_deg   Geodetic longitude (deg, East positive).
 * @param[in]  decimal_year    Time as decimal year (e.g., 2025.0).
 * @param[out] decl_deg        Declination D (deg). Set to NaN if H < 100 nT.
 * @param[out] incl_deg        Inclination I (deg, down-positive).
 * @param[out] F_nT            Total intensity F (nT).
 * @param[out] GV_deg          Grid variation (deg). -999.0 at low latitudes.
 * @param[out] X_nT            North component (nT) = F cos(D) cos(I).
 * @param[out] Y_nT            East  component (nT) = F sin(D) cos(I).
 * @param[out] Z_nT            Down  component (nT) = F sin(I).
 * @param[out] H_nT            Horizontal intensity (nT) = F cos(I).
 *
 * @return WMM_OK (0) or WMM_WARN(1) when horizontal field is weak.
 */
WMM_MAG_API int wmm_point(double altitude_km, double latitude_deg, double longitude_deg, double decimal_year,
              double *decl_deg, double *incl_deg, double *F_nT, double *GV_deg,
              double *X_nT, double *Y_nT, double *Z_nT, double *H_nT);

/**
 * @brief Compute elements and 1-year secular change at a point.
 *
 * @param[in]  altitude_km, latitude_deg, longitude_deg, decimal_year  Position/time.
 * @param[out] decl_deg, incl_deg, F_nT, GV_deg     Main elements.
 * @param[out] X_nT, Y_nT, Z_nT, H_nT               NED & horizontal.
 * @param[out] Xdot_nTyr, Ydot_nTyr, Zdot_nTyr      dX/dt, dY/dt, dZ/dt (nT/yr).
 * @param[out] Hdot_nTyr, Fdot_nTyr                 dH/dt, dF/dt (nT/yr).
 * @param[out] Ddot_degyr, Idot_degyr               dD/dt, dI/dt (deg/yr).
 *
 * @return WMM_OK (0) or WMM_WARN(1).
 *
 * @details
 * Evaluates values at t and at t+1 year internally, then differences them.
 */
WMM_MAG_API int wmm_point_change(double altitude_km, double latitude_deg, double longitude_deg, double decimal_year,
                     double *decl_deg, double *incl_deg, double *F_nT, double *GV_deg,
                     double *X_nT, double *Y_nT, double *Z_nT, double *H_nT,
                     double *Xdot_nTyr, double *Ydot_nTyr, double *Zdot_nTyr,
                     double *Hdot_nTyr, double *Fdot_nTyr,
                     double *Ddot_degyr, double *Idot_degyr);

/*======================================================================*
 *  Utilities
 *======================================================================*/

/**
 * @brief Compute cartesian components from (D,I,F).
 *
 * @param[in]  decl_deg  Declination (deg).
 * @param[in]  incl_deg  Inclination (deg).
 * @param[in]  F_nT      Total intensity F (nT).
 * @param[out] X_nT      North (nT) = F cosD cosI.
 * @param[out] Y_nT      East  (nT) = F sinD cosI.
 * @param[out] Z_nT      Down  (nT) = F sinI.
 * @param[out] H_nT      Horizontal (nT) = F cosI.
 */
WMM_MAG_API void wmm_components(double decl_deg, double incl_deg, double F_nT,
                    double *X_nT, double *Y_nT, double *Z_nT, double *H_nT);

/**
 * @brief Provide representative 1-σ uncertainties (WMM-2025 guidance).
 *
 * @param[in]  H_nT     Horizontal intensity (nT).
 * @param[out] uF_nT    σ(F)  (nT)
 * @param[out] uH_nT    σ(H)  (nT)
 * @param[out] uX_nT    σ(X)  (nT)
 * @param[out] uY_nT    σ(Y)  (nT)
 * @param[out] uZ_nT    σ(Z)  (nT)
 * @param[out] uI_deg   σ(I)  (deg)
 * @param[out] uD_deg   σ(D)  (deg) ~ sqrt(Do^2 + (K/H)^2)
 *
 * @note At very low H (e.g., near magnetic poles), angle uncertainties
 *       grow large and may be reported as NaN.
 */
WMM_MAG_API void wmm_error(double H_nT,
               double *uF_nT, double *uH_nT,
               double *uX_nT, double *uY_nT, double *uZ_nT,
               double *uI_deg, double *uD_deg);

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* WMM_MAG_H */
