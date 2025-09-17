/*======================================================================*
 *  File        : wmm_mag.c
 *  Brief       : World Magnetic Model (WMM-2025) public API (C)
 *  Author      : Nova Dynamics (novadynamics.kr@gmail.com)
 *  License     : MIT
 *
 *  Notes
 *  - This code is rewritten based on the WMM2020 Legacy C code.
 *      * https://www.ngdc.noaa.gov/geomag/WMM/data/WMM2020/WMM2020LegacyC.zip
 *      * Author: NCEI Geomagnetism Team
 *      * Date: December 2024
 *  - Units:
 *      * latitude/longitude: degrees (geodetic; East/North positive)
 *      * altitude          : kilometers above WGS84 ellipsoid
 *      * time              : decimal year (e.g., 2025.0)
 *      * field outputs     : nT; angles in degrees
 *======================================================================*/

#include "wmm_mag.h"
#include "wmm_mag_struct.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @name WMM 2025 representative 1-σ uncertainties
 * @details
 * The following are globally representative 1-σ uncertainty figures aligned
 * with WMM‑2025 guidance (units matched below).
 * - Units:
 *   - F, H, X, Y, Z : nT
 *   - I (inclination), D (declination) : degrees
 * - Declination uncertainty model: \f$ \sigma_D \approx \sqrt{D_{off}^2 + (D_{coef}/H)^2} \f$,
 *   i.e., heading uncertainty grows as horizontal strength H decreases.
 * @{ */
#define WMM_UNCERTAINTY_F 138
#define WMM_UNCERTAINTY_H 133
#define WMM_UNCERTAINTY_X 137
#define WMM_UNCERTAINTY_Y 89
#define WMM_UNCERTAINTY_Z 141
#define WMM_UNCERTAINTY_I 0.20
#define WMM_UNCERTAINTY_D_OFFSET 0.26
#define WMM_UNCERTAINTY_D_COEF 5417
/** @} */


/*======================================================================*
 *  Constants / Return codes
 *======================================================================*/

/**
 * @def WMM_MAXORD_DEFAULT
 * @brief Default maximum degree/order of the WMM (12 for WMM-2025).
 */
#ifndef WMM_MAXORD_DEFAULT
#define WMM_MAXORD_DEFAULT 12
#endif

/**
 * @name Warning flags (return value from wmm_point / wmm_point_change)
 * @note Current implementation returns 0 or 1 (generic warning).
 *       If you extend it to bit flags, consider:
 *       - WMM_WARN_LOW_H(5000 nT), WMM_WARN_VERY_LOW_H(1000 nT), WMM_WARN_POLE(100 nT)
 * @{ */
#define WMM_OK          0  /**< No warning */
#define WMM_WARN        1  /**< Generic warning (low-H or pole conditions) */
/** @} */


/* ---------- Globals ---------- */
static WMMState gWMM;
static WMMVersion gWMMVer;

/**
 * @brief for improved portablility. IEEE: only NaN is not equal to itself
 * @param d double value
 * @return 1 if d is NaN, else 0
 */
static int my_isnan(double d)
{
    return (d != d);
}

/**
 * @brief Convert calendar date/time to decimal year (UTC).
 *
 * @param year          Year (e.g., 2027)
 * @param month         Month [1..12]. If 0, returns @p year unchanged (legacy behavior).
 * @param day           Day within month [1..MonthDays]
 * @param hour          Hour [0..23]
 * @param minute        Minute [0..59]
 * @param[out] decimal_year  Decimal year result (e.g., 2027.500…)
 *
 * @return 1 on success, 0 on invalid input.
 *
 * @note If @p month == 0, the function returns @p year and ignores
 *       @p day/@p hour/@p minute (kept for backward compatibility).
 *
 * @details
 * - Leap-year rules: Gregorian (divisible by 400, or divisible by 4 but not by 100).
 * - Year fraction (DOY = day of year) is computed as:
 *   @code
 *   ( (DOY - 1) + hour/24 + minute/1440 ) / (365 or 366)
 *   @endcode
 */
int dec_year(int year, int month, int day, int hour, int minute, double *decimal_year)
{
    int temp = 0; /* Day-of-year accumulator (1-based later) */
    int MonthDays[13];
    int ExtraDay = 0;
    int i;

    double fractional_day = 0.0;

    if (month == 0) {
        /* Legacy behavior: just the year, ignore day/time */
        if (decimal_year) *decimal_year = year;
        return 1;
    }

    /* Leap year? */
    if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) {
        ExtraDay = 1;
    }

    MonthDays[0]  = 0;
    MonthDays[1]  = 31;
    MonthDays[2]  = 28 + ExtraDay;
    MonthDays[3]  = 31;
    MonthDays[4]  = 30;
    MonthDays[5]  = 31;
    MonthDays[6]  = 30;
    MonthDays[7]  = 31;
    MonthDays[8]  = 31;
    MonthDays[9]  = 30;
    MonthDays[10] = 31;
    MonthDays[11] = 30;
    MonthDays[12] = 31;

    /* Basic validation */
    if (month < 1 || month > 12) return 0;
    if (day   < 1 || day   > MonthDays[month]) return 0;
    if (hour  < 0 || hour  > 23) return 0;
    if (minute< 0 || minute> 59) return 0;

    /* Sum days before this month */
    for (i = 1; i <= month; i++) {
        temp += MonthDays[i - 1];
    }
    /* Add current day (1-based DOY) */
    temp += day;

    /* Fraction within the day */
    fractional_day = (hour / 24.0) + (minute / 1440.0);

    /* Decimal year */
    if (decimal_year) {
        *decimal_year = year + (((temp - 1) + fractional_day) / (365.0 + ExtraDay));
    }

    return 1;
}


/**
 * @brief Reset WMM internal state.
 * @param[out] S     State pointer
 * @param[in]  maxdeg Maximum degree/order (typically 12)
 * @return true on success
 * @details
 * - Zero-initialize arrays, set WGS84/WMM constants, set radians factor,
 *   and invalidate last-evaluation caches.
 */
static bool wmm_reset(WMMState *S, int maxdeg)
{
    /* Sets default values for WMM subroutines */
    memset(S, 0, sizeof(*S));
    S->maxord = maxdeg;

    /* base trig, normalization seeds */
    S->sp[0] = 0.0;
    S->cp[0] = 1.0;
    S->snorm[0] = 1.0;
    S->pp[0] = 1.0;
    S->dp[0][0] = 0.0;

    /* WGS84 / model constants */
    S->a  = 6378.137;      /* semi-major axis (km) */
    S->b  = 6356.7523142;  /* semi-minor axis (km) */
    S->re = 6371.2;        /* reference radius (km) */

    S->a2 = S->a*S->a;     S->b2 = S->b*S->b;     S->c2 = S->a2 - S->b2;
    S->a4 = S->a2*S->a2;   S->b4 = S->b2*S->b2;   S->c4 = S->a4 - S->b4;
    S->pi = M_PI;
    S->dtr = S->pi / 180.0;

    S->c[0][0]  = 0.0;
    S->cd[0][0] = 0.0;

    S->otime = S->oalt = S->olat = S->olon = -1e9;
    return true;
}

/**
 * @brief Load Gauss coefficients and secular variation from a COF file.
 * @param[out] S        State struct
 * @param[in]  cof_path COF file path
 * @param[out] ver      Optional version metadata (nullable)
 * @return true on success
 * @details
 * - Line 1: "<epoch> <model-name>"
 * - Line 2: (optional) Release date "MM/DD/YYYY"
 * - Following: rows of n m g h dg dh (until sentinel 9999)
 * - After loading, you must call @ref wmm_init_state.
 */
bool wmm_load_cof(WMMState *S, const char *cof_path, WMMVersion *ver)
{
    /* --------------------------------------------------------------------
    * COF parsing & storage layout (WMM-2025)
    *
    * The NOAA/BGS WMM COF file provides:
    *    n  m   g_nm   h_nm   dg_nm   dh_nm
    * where g/h are main-field Gauss coefficients [nT] at model epoch
    * and dg/dh are their linear secular-variation terms [nT/yr].
    *
    * This loader stores:
    *   c[m][n]   <- g_nm  (cosine terms)
    *   c[n][m-1] <- h_nm  (sine  terms)   for m>0
    * and similarly for cd[.][.] (time derivatives),
    * matching the classic NOAA WMM reference implementation scheme.
    *
    * After loading, call wmm_init() to convert coefficients into
    * Schmidt semi-normalized space used by the synthesis.
    *
    * Reference:
    *  - NCEI WMM coefficients & tech report links (WMM2025).
    * ------------------------------------------------------------------ */

    FILE *f = fopen(cof_path, "r");
    if (ver) {
        memset(ver, 0, sizeof(*ver));
        ver->epoch = NAN;
        ver->model_name[0] = '\0';
        ver->release_date[0] = '\0';
    }
    if (f == NULL) {
        fprintf(stderr, "Error opening model file %s\n", cof_path);
        return false;
    }

    /* header line: "<epoch> <model-name>" */
    char line[256];
    char model[32] = "";
    if (!fgets(line, sizeof(line), f) ||
        sscanf(line, "%lf%31s", &S->epoch, model) < 1) {
        fprintf(stderr, "Invalid header in model file %s\n", cof_path);
        fclose(f);
        return false;
    }
    if (ver) {
        ver->epoch = S->epoch;
        if (model[0] != '\0') {
            strncpy(ver->model_name, model, sizeof(ver->model_name)-1);
            ver->model_name[sizeof(ver->model_name)-1] = '\0';
        }
    }

    /* Try read 2nd line as release date or first coeff line */
    bool prefetched = false;
    if (fgets(line, sizeof(line), f)) {
        int slash = 0;
        for (char *p = line; *p; ++p) slash += (*p == '/');
        int n, m; double g, h, dg, dh;
        bool looks_like_coeff =
            (sscanf(line, "%d%d%lf%lf%lf%lf", &n, &m, &g, &h, &dg, &dh) == 6);
        if (slash >= 2 && !looks_like_coeff) {
            if (ver) {
                size_t len = strcspn(line, "\r\n");
                if (len >= sizeof(ver->release_date)) len = sizeof(ver->release_date)-1;
                memcpy(ver->release_date, line, len);
                ver->release_date[len] = '\0';
            }
        } else {
            prefetched = true;
        }
    }

    /* coefficients until sentinel '9999' */
    for (;;) {
        if (!prefetched) {
            if (fgets(line, sizeof(line), f) == NULL) break;
        } else {
            prefetched = false;
        }

        char c_new[5] = {0};
        for (int i = 0; i < 4 && line[i] != '\0'; ++i) c_new[i] = line[i];
        if (strcmp("9999", c_new) == 0) break;

        int n, m; double gnm, hnm, dgnm, dhnm;
        if (sscanf(line, "%d%d%lf%lf%lf%lf",
                   &n, &m, &gnm, &hnm, &dgnm, &dhnm) != 6) {
            continue; /* skip headers/comments */
        }

        //  gnm : C - Gauss coefficients of main geomagnetic model (nT)
        //  hnm : C - Gauss coefficients of main geomagnetic model (nT)
        // dgnm : CD - Gauss coefficients of secular geomagnetic model (nT/yr)
        // dhnm : CD - Gauss coefficients of secular geomagnetic model (nT/yr)
        if (n <= S->maxord && m >= 0 && m <= n) {
            S->c[m][n]  = gnm;
            S->cd[m][n] = dgnm;
            if (m != 0) {
                S->c[n][m-1]  = hnm;
                S->cd[n][m-1] = dhnm;
            }
        }
    }

    fclose(f);
    return true;
}

void wmm_version(double *epoch,
                 char *model_name, size_t model_name_cap,
                 char *release_date, size_t release_date_cap)
{
    if (epoch) {
        *epoch = gWMMVer.epoch;
    }

    if (model_name && model_name_cap > 0) {
        if (gWMMVer.model_name[0] != '\0') {
            strncpy(model_name, gWMMVer.model_name, model_name_cap - 1);
            model_name[model_name_cap - 1] = '\0';
        } else {
            model_name[0] = '\0';
        }
    }

    if (release_date && release_date_cap > 0) {
        if (gWMMVer.release_date[0] != '\0') {
            strncpy(release_date, gWMMVer.release_date, release_date_cap - 1);
            release_date[release_date_cap - 1] = '\0';
        } else {
            release_date[0] = '\0';
        }
    }
}


/**
 * @brief Initialize state: Schmidt normalization and Legendre recursion helpers.
 * @param[in,out] S  State struct
 * @return true on success
 * @details
 * - Convert g/h (and dg/dh) to Schmidt semi-normalized form for internal use.
 * - Prepare associated Legendre recursion coefficients and helper tables.
 * - Invalidate last-evaluation caches.
 */
bool wmm_init_state(WMMState* S)
{
    /* --------------------------------------------------------------------
    * Schmidt semi-normalization & recursion tables
    *
    * WMM uses Schmidt "semi-normalized" associated Legendre functions.
    * Here we (1) pre-compute normalization factors snorm[..],
    * (2) scale g/h and their time-derivatives accordingly, and
    * (3) prepare recursion helpers k[ m ][ n ], fn[n]=n+1, fm[m]=m.
    *
    * This makes the runtime synthesis loop (in wmm_geomag) numerically
    * stable and consistent with the COF's intended normalization.
    *
    * Key ideas:
    *  - snorm indexes map (n,m) into a flat array: snorm[n + m*13]
    *  - for m>0, sine-terms are stored at c[n][m-1] (legacy layout)
    *  - k[m][n] holds the recursion coefficient used in P_n^m recurrences
    *
    * Background:
    *  - Potential V expanded to degree/order 12, Schmidt semi-normalized
    *    P̄_n^m( cosθ ).  B = -∇V is synthesized from Br, Bθ, Bφ.
    * ------------------------------------------------------------------ */

    double *p = S->snorm;
    p[0] = 1.0;
    S->fm[0] = 0.0;

    for (int n = 1; n <= S->maxord; n++) {
        int j = 2;
        p[n] = p[n-1] * (double)(2*n-1) / (double)n;
        for (int m = 0; m <= n; m++) {
            S->k[m][n] = (double)(((n-1)*(n-1)) - (m*m))
                       / (double)((2*n-1)*(2*n-3));
            if (m > 0) {
                double flnmj = (double)((n - m + 1) * j) / (double)(n + m);
                S->snorm[n + m*13] = S->snorm[n + (m-1)*13] * sqrt(flnmj);
                j = 1;
                /* h_{n,m} is stored at [n][m-1] */
                S->c[n][m-1]  = S->snorm[n + m*13] * S->c[n][m-1];
                S->cd[n][m-1] = S->snorm[n + m*13] * S->cd[n][m-1];
            }
            S->c[m][n]  = S->snorm[n + m*13] * S->c[m][n];
            S->cd[m][n] = S->snorm[n + m*13] * S->cd[m][n];
        }
        S->fn[n] = (double)(n + 1);
        S->fm[n] = (double)n;
    }
    S->k[1][1] = 0.0;

    S->otime = S->oalt = S->olat = S->olon = -1e9;
    return true;
}

/**
 * @brief Compute geomagnetic elements (Decl, Incl, F, GV) at a location/time.
 * @param[in,out] S   State
 * @param[in]     alt Altitude (km, WGS84 ellipsoid)
 * @param[in]     glat Geodetic latitude (deg)
 * @param[in]     glon Geodetic longitude (deg, East positive)
 * @param[in]     time Decimal year (dt = time - epoch)
 * @param[out]    dec  Declination D (deg, east-positive)
 * @param[out]    dip  Inclination I (deg, down-positive)
 * @param[out]    f    Total intensity F (nT)
 * @param[out]    gv   Grid variation (deg; meaningful at |lat|≥55°, otherwise -999.0)
 * @details
 * - Time-update g/h → geodetic→geocentric conversion → Schmidt-normalized Legendre
 *   and derivative → spherical-harmonic synthesis → geocentric→geodetic local rotation.
 * - Cartesian components (X,Y,Z,H) are not returned here; use @ref wmm_components if needed.
 */
void wmm_geomag(WMMState* S,
                     double alt, double glat, double glon, double time,
                     double* dec, double* dip, double* f, double* gv)
{
    double dt = time - S->epoch;

    double rlon = glon * S->dtr;
    double rlat = glat * S->dtr;
    double srlon = sin(rlon), crlon = cos(rlon);
    double srlat = sin(rlat), crlat = cos(rlat);
    double srlat2 = srlat*srlat, crlat2 = crlat*crlat;

    S->sp[1] = srlon;
    S->cp[1] = crlon;

    /* === Geodetic (lat,lon,alt) -> Geocentric spherical (r, theta, lambda) ===
    *
    * Inputs:
    *   glat, glon in degrees (geodetic latitude/longitude, WGS84)
    *   alt in km above WGS84 ellipsoid
    *
    * We convert geodetic latitude to geocentric colatitude theta via the
    * standard WGS84 geometry. Symbols:
    *   r   : geocentric radius (km)
    *   ct  : cos(theta) = cos(geocentric colatitude)
    *   st  : sin(theta)
    *   ca  : cos(geodetic->geocentric tilt)
    *   sa  : sin(geodetic->geocentric tilt)
    *
    * Notes:
    *   - WMM reference radius Re=6371.2 km is used in synthesis, not WGS84 a.
    *   - 'ca' and 'sa' later rotate geocentric (Br,Bθ,Bφ) into geodetic local
    *     North-East-Down components (X=N, Y=E, Z=Down).
    *
    * Ref. concepts: WMM-2025 report (NCEI), potential expansion in spherical
    * harmonics with Re=6371.2 km and Schmidt normalization.
    * -------------------------------------------------------------------- */
    double ct = S->ct;
    double st = S->st;
    double r = S->r;
    double ca = S->ca;
    double sa = S->sa;
    if (alt != S->oalt || glat != S->olat) {
        double q  = sqrt(S->a2 - S->c2 * srlat2);
        double q1 = alt * q;
        double q2 = ((q1 + S->a2) / (q1 + S->b2));
        q2 = q2*q2;
        ct = srlat / sqrt(q2 * crlat2 + srlat2);
        st = sqrt(fmax(0.0, 1.0 - ct*ct));
        double r2 = (alt*alt) + 2.0*q1 + (S->a4 - S->c4 * srlat2) / (q*q);
        r  = sqrt(r2);
        double d  = sqrt(S->a2 * crlat2 + S->b2 * srlat2);
        ca = (alt + d) / r;
        sa = (S->c2 * crlat * srlat) / (r * d);

        S->ct = ct; S->st = st; S->r = r; S->ca = ca; S->sa = sa;
    }

    if (glon != S->olon) {
        for (int m = 2; m <= S->maxord; m++) {
            S->sp[m] = S->sp[1] * S->cp[m-1] + S->cp[1] * S->sp[m-1];
            S->cp[m] = S->cp[1] * S->cp[m-1] - S->sp[1] * S->sp[m-1];
        }
    }

    /* === Schmidt-normalized associated Legendre P̄_n^m and dP̄/dθ (on the fly) ===
    *
    * We compute P̄ and their derivatives via recursions driven by:
    *   snorm[n + m*13]  : accumulated Schmidt normalization
    *   k[m][n]          : recurrence coefficient
    *   dp[m][n]         : stores ∂P̄_n^m / ∂θ (theta-derivative)
    *
    * Cases handled:
    *   - Diagonal (n==m):   P̄_n^n  from P̄_{n-1}^{n-1} via sin(theta)
    *   - Subdiagonal:       P̄_n^{n-1} from P̄_{n-1}^{n-1} via cos(theta)
    *   - General m<=n-2:    full two-step recurrence using k[m][n]
    *
    * Numerics:
    *   - We avoid evaluating at exact st=0 by deferring bp (Bφ) handling.
    *   - dp stores d/dθ of already Schmidt-normalized values.
    * -------------------------------------------------------------------- */

    /* === Spherical harmonic synthesis (main field only; no external terms) ===
    *
    * Magnetic potential:
    *   V(r,θ,λ) = Re * Σ_{n=1..N} (Re/r)^{n+1} Σ_{m=0..n}
    *              [ ḡ_n^m cos(mλ) + h̄_n^m sin(mλ) ] P̄_n^m( cosθ )
    * Field:
    *   B = -∇V  ⇒ components in spherical basis (r, θ, φ):
    *     Br     = + Σ (n+1) (Re/r)^{n+2} T * P̄_n^m
    *     Bθ     = - Σ        (Re/r)^{n+2} T * ∂P̄_n^m/∂θ
    *     Bφ     = + Σ        (Re/r)^{n+2} (m/ sinθ) * Tφ * P̄_n^m
    *   where T   = ḡ_n^m cos(mλ) + h̄_n^m sin(mλ)
    *         Tφ  = ḡ_n^m sin(mλ) - h̄_n^m cos(mλ)  (longitude derivative)
    *
    * Implementation notes:
    *   - aor=Re/r; ar accumulates (Re/r)^{n+2}
    *   - cml/sml are cos(mλ)/sin(mλ) via cheap recurrences
    *   - Polar case (st≈0): use bpp accumulator then set bp=bpp (avoids /0)
    *
    * Time update (secular variation):
    *   ḡ_n^m(t) = ḡ_n^m(epoch) + dt * ḡ̄_n^m, similarly for h̄
    *   (dt in years; WMM assumes linear SV over 5-year validity window)
    * -------------------------------------------------------------------- */

    double aor = S->re / r; /* 1/r scale */
    double ar  = aor * aor;
    double br = 0.0, bt = 0.0, bp = 0.0, bpp = 0.0;

    for (int n = 1; n <= S->maxord; n++) {
        ar *= aor;
        for (int m = 0; m <= n; m++) {

            if (alt != S->oalt || glat != S->olat) {
                if (n == m) {
                    S->snorm[n + m*13] = st * S->snorm[n-1 + (m-1)*13];
                    S->dp[m][n] = st * S->dp[m-1][n-1] + ct * S->snorm[n-1 + (m-1)*13];
                } else if (n == 1 && m == 0) {
                    S->snorm[n + m*13] = ct * S->snorm[n-1 + m*13];
                    S->dp[m][n] = ct * S->dp[m][n-1] - st * S->snorm[n-1 + m*13];
                } else if (n > 1 && n != m) {
                    if (m > n-2) { S->snorm[n-2 + m*13] = 0.0; S->dp[m][n-2] = 0.0; }
                    S->snorm[n + m*13] = ct * S->snorm[n-1 + m*13] - S->k[m][n] * S->snorm[n-2 + m*13];
                    S->dp[m][n] = ct * S->dp[m][n-1] - st * S->snorm[n-1 + m*13] - S->k[m][n] * S->dp[m][n-2];
                }
            }

            /* time-adjusted Gauss coeffs */
            if (time != S->otime) {
                S->tc[m][n] = S->c[m][n] + dt * S->cd[m][n];
                if (m != 0) S->tc[n][m-1] = S->c[n][m-1] + dt * S->cd[n][m-1];
            }

            double par = ar * S->snorm[n + m*13];
            double temp1, temp2;
            if (m == 0) {
                temp1 = S->tc[m][n] * S->cp[m];
                temp2 = S->tc[m][n] * S->sp[m];
            } else {
                temp1 = S->tc[m][n] * S->cp[m] + S->tc[n][m-1] * S->sp[m];
                temp2 = S->tc[m][n] * S->sp[m] - S->tc[n][m-1] * S->cp[m];
            }

            bt = bt - ar * temp1 * S->dp[m][n];
            bp += (S->fm[m] * temp2 * par);
            br += (S->fn[n] * temp1 * par);

            /* polar special case */
            if (st == 0.0 && m == 1) {
                if (n == 1) S->pp[n] = S->pp[n-1];
                else        S->pp[n] = ct * S->pp[n-1] - S->k[m][n] * S->pp[n-2];
                double parp = ar * S->pp[n];
                bpp += (S->fm[m] * temp2 * parp);
            }
        }
    }

    if (fabs(st) < 1e-12) bp = bpp;
    else                  bp /= st;

    /* === Rotate (Br,Bθ,Bφ) to local geodetic NED and compute elements ====
    *
    * Spherical basis (geocentric):
    *   r: outward, θ: southward, φ: eastward
    * Map to local geocentric (North, East, Down):
    *   N_gc = -Bθ,  E_gc =  Bφ,  D_gc = -Br
    * Then tilt from geocentric to geodetic by (ca,sa):
    *   X (North) =  N_gc*ca + D_gc*sa
    *   Y (East)  =  E_gc
    *   Z (Down)  = -N_gc*sa + D_gc*ca
    *
    * Field elements:
    *   H = sqrt(X^2 + Y^2),  F = sqrt(H^2 + Z^2)
    *   D = atan2(Y, X) [deg, east-positive],  I = atan2(Z, H) [deg, down+]
    *
    * Grid Variation (GV):
    *   Only meaningful at |lat|≥55°, derived from D and longitude quadrant.
    * -------------------------------------------------------------------- */


    /* rotate to geodetic components */
    double bx = -bt * ca - br * sa;
    double by =  bp;
    double bz =  bt * sa - br * ca;

    double bh = hypot(bx, by);
    *f  = hypot(bh, bz);
    *dec = atan2(by, bx) / S->dtr;
    *dip = atan2(bz, bh) / S->dtr;

    /* grid variation (high latitudes) */
    *gv = -999.0;
    if (fabs(glat) >= 55.0) {
        if (glat > 0.0 && glon >= 0.0) *gv = *dec - glon;
        if (glat > 0.0 && glon <  0.0) *gv = *dec + fabs(glon);
        if (glat < 0.0 && glon >= 0.0) *gv = *dec + glon;
        if (glat < 0.0 && glon <  0.0) *gv = *dec - fabs(glon);
        if (*gv >  180.0) *gv -= 360.0;
        if (*gv < -180.0) *gv += 360.0;
    }

    S->otime = time;
    S->oalt  = alt;
    S->olat  = glat;
    S->olon  = glon;
}

/**
 * @brief (Wrapper) Initialize global WMM state: Reset → Load COF → Init.
 * @param[in] cof    COF file path
 */
void wmm_init(const char *cof) {
    int maxdeg = 12;
    if (!wmm_reset(&gWMM, maxdeg)) {
        fprintf(stderr, "WMM reset failed\n");
        exit(EXIT_FAILURE);
    }
    if (!wmm_load_cof(&gWMM, cof, &gWMMVer)) {
        fprintf(stderr, "WMM load failed for %s\n", cof);
        exit(EXIT_FAILURE);
    }
    if (!wmm_init_state(&gWMM)) {
        fprintf(stderr, "WMM init failed\n");
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Compute Cartesian components (X,Y,Z,H) from (D, I, F).
 * @param[in]  dec_deg Declination D (deg)
 * @param[in]  dip_deg Inclination I (deg)
 * @param[in]  f       Total intensity F (nT)
 * @param[out] x       X (North, nT) = F cos(D) cos(I)
 * @param[out] y       Y (East, nT)  = F sin(D) cos(I)
 * @param[out] z       Z (Down, nT)  = F sin(I)
 * @param[out] h       H (Horizontal, nT) = F cos(I)
 */
void wmm_components(double dec_deg, double dip_deg, double f,
                           double *x, double *y, double *z, double *h)
{
    const double rTd = M_PI / 180.0;
    //const double rTd=0.017453292;
    const double cdec = cos(dec_deg * rTd);
    const double sdec = sin(dec_deg * rTd);
    const double cdip = cos(dip_deg * rTd);
    const double sdip = sin(dip_deg * rTd);

    if (x) *x = f * (cdec * cdip);
    if (y) *y = f * (cdip * sdec);
    if (z) *z = f *  sdip;
    if (h) *h = f *  cdip;
}

/**
 * @brief Compute secular change by differencing two epochs separated by 1 year.
 * @param[in]  f1     Baseline total F1 at t (nT)
 * @param[in]  x1,y1,z1,h1  Baseline components at t (nT)
 * @param[in]  x2,y2,z2     Components at t+1 (nT)
 * @param[out] xdot, ydot, zdot (nT/yr)
 * @param[out] hdot, fdot        (nT/yr)
 * @param[out] ddot, idot        (deg/yr) — annual change of D and I
 */
void wmm_changes(double f1, double x1, double y1, double z1, double h1,
                        double x2, double y2, double z2,
                        double *xdot, double *ydot, double *zdot, double *hdot, double *fdot,
                        double *ddot, double *idot)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;

    double dh = (x1 * dx + y1 * dy) / h1;
    double df = (x1 * dx + y1 * dy + z1 * dz) / f1;
    *ddot = 180.0 / M_PI * (x1 * dy - y1 * dx) / (h1 * h1);
    *idot = 180.0 / M_PI * (h1 * dz - z1 * dh) / (f1 * f1);

    *xdot  = dx;
    *ydot  = dy;
    *zdot  = dz;

    *hdot  = dh;
    *fdot = df;
}

/**
 * @brief Estimate 1‑year change by evaluating at (t+1 yr) and differencing.
 * @param[in]  S                    State
 * @param[in]  alt,lat,lon,time     Location/time
 * @param[in]  f,x,y,z,h            Baseline total/components at time t
 * @param[out] xdot,ydot,zdot       Component rates (nT/yr)
 * @param[out] hdot,fdot            H and F rates (nT/yr)
 * @param[out] ddot,idot            D and I rates (deg/yr)
 * @details Internally re-evaluates at t+1 year, then calls @ref wmm_changes.
 *          Mirrors the WMM linear-SV assumption over a 5‑year validity window.
 */
void wmm_annual_change(WMMState *S,
                              double alt, double lat, double lon, double time,
                              double f, double x, double y, double z, double h,
                              double *xdot, double *ydot, double *zdot, double *hdot,
                              double *fdot, double *ddot, double *idot)
{
    // value at t+1 year
    double dec2, dip2, f2, gv2;
    wmm_geomag(S, alt, lat, lon, time + 1.0, &dec2, &dip2, &f2, &gv2);

    double x2, y2, z2, h2;
    wmm_components(dec2, dip2, f2, &x2, &y2, &z2, &h2);

    wmm_changes(f, x, y, z, h,
                x2, y2, z2,
                xdot, ydot, zdot, hdot, fdot,
                ddot, idot);
}

/**
 * @brief Provide representative 1‑σ uncertainties (D depends on H).
 * @param[in]  h      Horizontal component H (nT)
 * @param[out] uf     σ(F)  (nT)
 * @param[out] uh     σ(H)  (nT)
 * @param[out] ux,uy,uz  σ(X), σ(Y), σ(Z) (nT)
 * @param[out] uincl  σ(I)  (deg)
 * @param[out] udecl  σ(D)  (deg) — \f$ \sqrt{D_{off}^2 + (D_{coef}/H)^2} \f$
 * @note If H is very small or NaN, angle uncertainties are set to NaN.
 * @details Returns global-average representative 1‑σ values for quick guidance.
 */
void wmm_error(double h,
               double *uf, double *uh,
               double *ux, double *uy, double *uz,
               double *uincl, double *udecl)
{
    double decl_variable;
    double decl_constant;

    if (my_isnan(h) || (h < 100.0))
    {
        /* Magnetic pole (very small H): heading/angle components unreliable */
        if (uf)    *uf    = NAN;
        if (uh)    *uh    = NAN;
        if (ux)    *ux    = NAN;
        if (uy)    *uy    = NAN;
        if (uz)    *uz    = NAN;
        if (uincl) *uincl = NAN;
        if (udecl) *udecl = NAN;
        return;
    }

    if (uf)    *uf    = WMM_UNCERTAINTY_F;
    if (uh)    *uh    = WMM_UNCERTAINTY_H;
    if (ux)    *ux    = WMM_UNCERTAINTY_X;
    if (uz)    *uz    = WMM_UNCERTAINTY_Z;
    if (uincl) *uincl = WMM_UNCERTAINTY_I;
    if (uy)    *uy    = WMM_UNCERTAINTY_Y;

    decl_variable = (WMM_UNCERTAINTY_D_COEF / h);
    decl_constant = (WMM_UNCERTAINTY_D_OFFSET);
    if (udecl) {
        *udecl = sqrt(decl_constant*decl_constant + decl_variable*decl_variable);
        if (*udecl > 180) *udecl = 180;
    }
}


/**
 * @brief Compute WMM elements and Cartesian components, and return a warning flag.
 *
 * Given position (lat/lon/alt) and time (decimal year), computes
 * - Magnetic elements: Declination (D), Inclination (I), Total intensity (F),
 *   and Grid Variation (GV).
 * - Cartesian components: X (North), Y (East), Z (Down), and H (Horizontal).
 * A warning is returned when H is weak (e.g., near the magnetic poles),
 * indicating heading uncertainty.
 *
 * @pre Global WMM state (gWMM) must be initialized (e.g., via wmm_init()).
 * @see wmm_geomag, wmm_components
 *
 * @param[in]  altitude     Altitude (km, WGS84 ellipsoid)
 * @param[in]  latitude     Latitude (deg, geodetic)
 * @param[in]  longitude    Longitude (deg, geodetic, East-positive)
 * @param[in]  decimal_year Decimal year (dt = decimal_year - epoch)
 * @param[out] decl         D (deg). Set to NaN if H < 100 nT
 * @param[out] incl         I (deg)
 * @param[out] f            F (nT)
 * @param[out] gv           GV (deg; -999.0 at low latitudes)
 * @param[out] x,y,z,h      X,Y,Z,H components (nT)
 *
 * @return Warning flag
 * @retval 0 No warning
 * @retval 1 Warning (any of the following):
 *         - H < 5000 nT : heading uncertainty increases
 *         - H < 1000 nT : heading uncertainty is very large
 *         - H < 100  nT : near pole; declination set to NaN
 *
 * @note Only a single-bit warning is returned; for tiered handling,
 *       inspect H at the call site.
 */
int wmm_point(double altitude, double latitude, double longitude, double decimal_year,
                      double *decl, double *incl, double *f, double *gv,
                      double *x, double *y, double *z, double *h)
{
    int warning = 0;

    wmm_geomag(&gWMM, altitude, latitude, longitude, decimal_year, decl, incl, f, gv); // compute elements
    wmm_components(*decl, *incl, *f, x, y, z, h); // compute X, Y, Z, H

    double h_val = *h;

    // Set warning flags for weak horizontal field
    if (h_val < 1000.0) warning = 1;
    else if (h_val < 5000.0) warning = 1;

    /* Pole handling (very weak H): declination undefined */
    if (h_val < 100.0) {
        *decl = NAN; // declination not defined
        warning = 1;
    }

    return warning;
}

/**
 * @brief Compute elements/components and annual change; return a warning flag.
 *
 * For the given position and time, computes D, I, F, GV and X,Y,Z,H, and then
 * uses a t→t+1 comparison to estimate annual changes of each quantity.
 *
 * @pre Global WMM state (gWMM) must be initialized.
 * @see wmm_point, wmm_annual_change
 *
 * @param[in]  altitude,latitude,longitude,decimal_year  Location/time
 * @param[out] decl,incl,f,gv   D,I,F,GV (deg,deg,nT,deg)
 * @param[out] x,y,z,h          X,Y,Z,H (nT)
 * @param[out] xdot,ydot,zdot   X,Y,Z rates (nT/yr)
 * @param[out] hdot,fdot        H,F rates (nT/yr)
 * @param[out] ddot,idot        D,I rates (deg/yr)
 * @return Warning flag (0=OK, 1=Warn)
 * @details Internally calls @ref wmm_point and then @ref wmm_annual_change.
 */
int wmm_point_change(double altitude, double latitude, double longitude, double decimal_year,
        double *decl, double *incl, double *f, double *gv,
        double *x, double *y, double *z, double *h,
        double *xdot, double *ydot, double *zdot, double *hdot, double *fdot,
        double *ddot, double *idot)
{
    int warning = 0;

    // Compute baseline at time t
    warning = wmm_point(altitude, latitude, longitude, decimal_year,
        decl, incl, f, gv, x, y, z, h);

    // Compute annual change via t+1 difference
    wmm_annual_change(&gWMM, altitude, latitude, longitude, decimal_year,
                      *f, *x, *y, *z, *h,
                      xdot, ydot, zdot, hdot, fdot, ddot, idot);

    return warning;
}
