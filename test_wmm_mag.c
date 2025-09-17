#include "wmm_mag.h"
#include "wmm_mag_struct.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * @brief Print calculation results in table format.
 * @param[in] dlat,dlon,alt,time Location/time
 * @param[in] mag        Output values (D, I, F, GV, X, Y, Z, H, and annual changes)
 * @param[in] mag_error  Uncertainty (1-σ) container
 *                       - In this implementation:
 *                         mag_error->X,Y,Z,H,F  : 1-σ for each component/total (nT)
 *                         mag_error->Incldot    : 1-σ for I (deg)
 *                         mag_error->Decldot    : 1-σ for D (deg)
 * @param[in] warn_H, warn_H_val               Warning flag/value for H<5000 nT
 * @param[in] warn_H_strong, warn_H_strong_val Warning flag/value for H<1000 nT
 * @param[in] warn_P                           Pole warning (if needed)
 */
void wmm_print_results(double dlat, double dlon, double alt, double time,
                              const MAGtype_GeoMagneticElements* mag,
                              const MAGtype_GeoMagneticElements* mag_error,
                              int warn_H, double warn_H_val,
                              int warn_H_strong, double warn_H_strong_val,
                              int warn_P)
{
    printf("\n Results For \n");
    printf("\n      Latitude: %7.2lf%c", fabs(dlat), (dlat<0)?'S':'N');
    printf("\n     Longitude: %7.2lf%c", fabs(dlon), (dlon<0)?'W':'E');
    printf("\n      Altitude: %8.2lf KM ABOVE WGS84 ELLIPSOID", alt);
    printf("\n Decimal Year : %6.1lf\n", time);

    printf("\n          Main Field    \t\t Secular Change   \t Uncertainty (+/-)");

    if (isnan(mag->X)) printf("\n X      =    NaN         \t   dX  = NaN");
    else               printf("\n X      =    %-+9.1lf nT\t   dX  = %-+8.1lf nT/yr   \t %-8.1lf nT",
                              mag->X, mag->Xdot, mag_error->X);

    if (isnan(mag->Y)) printf("\n Y      =    NaN         \t   dY  = NaN");
    else               printf("\n Y      =    %-+9.1lf nT\t   dY  = %-+8.1lf nT/yr   \t %-8.1lf nT",
                              mag->Y, mag->Ydot, mag_error->Y);

    printf("\n Z      =    %-+9.1lf nT\t   dZ  = %-+8.1lf nT/yr   \t %-8.1lf nT",
           mag->Z, mag->Zdot, mag_error->Z);

    if (isnan(mag->H)) printf("\n H      =    NaN         \t   dH  = NaN");
    else               printf("\n H      =    %-9.1lf nT\t   dH  = %-+8.1lf nT/yr   \t %-8.1lf nT",
                              mag->H, mag->Hdot, mag_error->H);

    printf("\n F      =    %-9.1lf nT\t   dF  = %-+8.1lf nT/yr   \t %-8.1lf nT",
           mag->F, mag->Fdot, mag_error->F);

    printf("\n I      =    %-+6.02lf    deg\t   dI  = %-+8.2lf deg/yr   \t %-8.2lf deg",
           mag->Incl, mag->Incldot, mag_error->Incldot);

    if (isnan(mag->Decl))
        printf("\n D      =    NaN         \t   dD  = NaN");
    else
        printf("\n D      =    %-+6.02lf    deg\t   dD  = %-8.2lf deg/yr   \t %-8.2lf deg",
               mag->Decl, mag->Decldot, mag_error->Decldot);

    if (warn_H) {
        printf("\n\nWarning: The horizontal field strength at this location is only %6.1lf nT\n", warn_H_val);
        printf("         Compass readings have large uncertainties where H < 5000 nT\n");
    }
    if (warn_H_strong) {
        printf("\n\nWarning: The horizontal field strength at this location is only %6.1lf nT\n", warn_H_strong_val);
        printf("         Compass readings have VERY LARGE uncertainties where H < 1000 nT\n");
    }
    if (warn_P) {
        printf("\n\nWarning: Location is at geographic pole where X, Y, and Decl are undefined\n");
    }

    printf("\n");
}


/** @brief One test vector (one line) */
typedef struct {
    int     lineno;   /**< @brief Original line number in input file (for debugging) */
    double  date;     /**< @brief Field 1: Decimal year */
    double  alt;      /**< @brief Field 2: Height above WGS84 ellipsoid (km) */
    double  lat;      /**< @brief Field 3: Geodetic latitude (deg) */
    double  lon;      /**< @brief Field 4: Geodetic longitude (deg) */
    double  x;        /**< @brief Field 5: X (nT) */
    double  y;        /**< @brief Field 6: Y (nT) */
    double  z;        /**< @brief Field 7: Z (nT) */
    double  h;        /**< @brief Field 8: H (nT) */
    double  f;        /**< @brief Field 9: F (nT) */
    double  incl;     /**< @brief Field10: Inclination (deg) */
    double  decl;     /**< @brief Field11: Declination (deg) */
    double  gv;       /**< @brief Field12: Grid Variation (deg) (can be NaN) */
    double  xdot;     /**< @brief Field13: Xdot (nT/yr) */
    double  ydot;     /**< @brief Field14: Ydot (nT/yr) */
    double  zdot;     /**< @brief Field15: Zdot (nT/yr) */
    double  hdot;     /**< @brief Field16: Hdot (nT/yr) */
    double  fdot;     /**< @brief Field17: Fdot (nT/yr) */
    double  idot;     /**< @brief Field18: Idot (deg/yr) */
    double  ddot;     /**< @brief Field19: Ddot (deg/yr) */
} WMMTestCase;

/**
 * @brief Parse a whitespace-separated token as double (accepts "NaN").
 * @param[in]  s  Token string
 * @param[out] v  Converted value
 * @return 1 on success, 0 on failure
 */
int parse_double_token(const char* s, double* v) {
    if (!s || !*s) return 0;
    /* strtod supports "nan"/"NaN" (C99); generally OK across implementations */
    char* endp = NULL;
    double val = strtod(s, &endp);
    if (endp == s) return 0; /* conversion failed */
    *v = val;
    return 1;
}

#include <ctype.h>
#include <string.h>

#if defined(_MSC_VER)
#define strtok_r(s, d, p) strtok_s(s, d, p)
#endif

#define NF 19
/**
 * @brief Parse one line into a 19-field test case.
 * @param[in]  line   Input line (may include newline)
 * @param[in]  lineno Original line number (for debugging)
 * @param[out] tc     Output struct
 * @return 1 on success, 0 on skip (comment/blank) or failure
 */
int parse_test_line(const char* line, int lineno, WMMTestCase* tc) {
    /* Skip if the line is blank or the first non-space character is '#' */
    const char* p = line;
    while (*p && isspace((unsigned char)*p)) ++p;
    if (*p == '\0' || *p == '#') return 0; /* skip */

    /* Tokenization */
    char buf[512];
    strncpy(buf, line, sizeof(buf)-1);
    buf[sizeof(buf)-1] = '\0';

    //const int NF = 19;
    const char* tok[NF];
    int nt = 0;
    char* saveptr = NULL;
    for (char* t = strtok_r(buf, " \t\r\n", &saveptr); t; t = strtok_r(NULL, " \t\r\n", &saveptr)) {
        if (*t == '#') break; /* trailing comment */
        if (nt < NF) tok[nt++] = t;
    }
    if (nt < NF) return 0; /* malformed */

    double vals[NF];
    for (int i = 0; i < NF; ++i) {
        if (!parse_double_token(tok[i], &vals[i])) return 0;
    }

    /* Mapping */
    memset(tc, 0, sizeof(*tc));
    tc->lineno = lineno;
    tc->date   = vals[0];
    tc->alt    = vals[1];
    tc->lat    = vals[2];
    tc->lon    = vals[3];
    tc->x      = vals[4];
    tc->y      = vals[5];
    tc->z      = vals[6];
    tc->h      = vals[7];
    tc->f      = vals[8];
    tc->incl   = vals[9];
    tc->decl   = vals[10];
    tc->gv     = vals[11];
    tc->xdot   = vals[12];
    tc->ydot   = vals[13];
    tc->zdot   = vals[14];
    tc->hdot   = vals[15];
    tc->fdot   = vals[16];
    tc->idot   = vals[17];
    tc->ddot   = vals[18];
    return 1;
}

/**
 * @brief Read test file and build array of cases.
 * @param[in]  path    File path (e.g., "WMM2025_TEST_VALUES.txt")
 * @param[out] cases   Dynamically allocated array (must free() on success)
 * @param[out] count   Number of cases
 * @return 1 on success, 0 on failure
 */
int load_wmm_test_values(const char* path, WMMTestCase** cases, size_t* count) {
    FILE* f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Failed to open %s\n", path);
        return 0;
    }
    size_t cap = 32, n = 0;
    WMMTestCase* arr = (WMMTestCase*)malloc(cap * sizeof(WMMTestCase));
    if (!arr) { fclose(f); return 0; }

    char line[512];
    int lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        ++lineno;
        WMMTestCase tc;
        if (!parse_test_line(line, lineno, &tc)) continue;
        if (n == cap) {
            cap *= 2;
            WMMTestCase* tmp = (WMMTestCase*)realloc(arr, cap * sizeof(WMMTestCase));
            if (!tmp) { free(arr); fclose(f); return 0; }
            arr = tmp;
        }
        arr[n++] = tc;
    }
    fclose(f);
    *cases = arr;
    *count = n;
    return 1;
}

/** @brief |a-b|<=tol check. NaN is treated specially (both NaN = "equal"). */
static int nearly_equal(double a, double b, double tol) {
    if (isnan(a) && isnan(b)) return 1;
    if (isnan(a) || isnan(b)) return 0;
    return fabs(a - b) <= tol;
}

/**
 * @brief Compare wmm_point_change() results to expected values for one case.
 * @param[in]  tc     Test vector
 * @param[in]  tol    Allowed tolerance (e.g., 0.1)
 * @param[out] ok     Overall pass/fail (1/0)
 * @details
 * - Compared fields: X,Y,Z,H,F, I, D, (GV*), Xdot,Ydot,Zdot,Hdot,Fdot, Idot, Ddot
 * - *GV is skipped if expected value is NaN or computed value is -999.0.
 */
static void check_one_case(const WMMTestCase* tc, double tol, int* ok) {
    double decl, incl, f, gv, x, y, z, h;
    double ax, ay, az, ah, af, adec, ainc;

    (void)wmm_point_change(tc->alt, tc->lat, tc->lon, tc->date,
            &decl, &incl, &f, &gv,
            &x, &y, &z, &h,
            &ax, &ay, &az, &ah, &af,
            &adec, &ainc);

    int pass = 1;

    /* Scalar comparison helper (prints a message on failure) */
    #define CHECK(lbl, got, exp) do { \
        if (!nearly_equal((got), (exp), tol)) { \
            fprintf(stderr, "  FAIL %-4s : got=% .6f  exp=% .6f  |Δ|=%.6f\n", \
                    (lbl), (got), (exp), fabs((got)-(exp))); \
            pass = 0; \
        } \
    } while(0)

    CHECK("X",     x,     tc->x);
    CHECK("Y",     y,     tc->y);
    CHECK("Z",     z,     tc->z);
    CHECK("H",     h,     tc->h);
    CHECK("F",     f,     tc->f);
    CHECK("I",     incl,  tc->incl);
    CHECK("D",     decl,  tc->decl);

    /* GV: skip comparison if expected is NaN or computed value is -999.0 (low latitude) */
    if (!(isnan(tc->gv) || gv == -999.0)) {
        CHECK("GV", gv, tc->gv);
    }

    CHECK("Xdot",  ax,    tc->xdot);
    CHECK("Ydot",  ay,    tc->ydot);
    CHECK("Zdot",  az,    tc->zdot);

    CHECK("Hdot",  ah,    tc->hdot);
    CHECK("Fdot",  af,    tc->fdot);

    CHECK("Idot",  ainc,  tc->idot);
    CHECK("Ddot",  adec,  tc->ddot);

    #undef CHECK
    *ok = pass;
}

/**
 * @brief Load test file and run/validate all cases.
 * @param[in] path  Test file path (e.g., "WMM2025_TEST_VALUES.txt")
 * @param[in] tol   Allowed tolerance (e.g., 0.1)
 * @return 0(success) / non-0(failure)
 * @note Global WMM state must be initialized (call wmm_init first).
 */
static int run_wmm_tests(const char* path, double tol) {
    WMMTestCase* cases = NULL;
    size_t n = 0;
    if (!load_wmm_test_values(path, &cases, &n)) {
        fprintf(stderr, "Could not load test cases from %s\n", path);
        return 1;
    }
    if (n == 0) {
        fprintf(stderr, "No test cases found in %s\n", path);
        free(cases);
        return 1;
    }

    size_t n_pass = 0;
    for (size_t i = 0; i < n; ++i) {
        const WMMTestCase* tc = &cases[i];
        int ok = 0;
        printf("\n[Case %zu] date=%.3f alt=%.1f km lat=%.2f lon=%.2f\n",
               i+1, tc->date, tc->alt, tc->lat, tc->lon);
        check_one_case(tc, tol, &ok);
        if (ok) {
            printf("  -> PASS (tol=%.3f)\n", tol);
            ++n_pass;
        } else {
            printf("  -> FAIL (tol=%.3f)  [line %d]\n", tol, tc->lineno);
        }
    }

    printf("\n== SUMMARY ==  %zu / %zu passed (tol=%.3f)\n", n_pass, n, tol);
    free(cases);
    return (n_pass == n) ? 0 : 2;
}

/**
 * @brief Sample main: compute and print WMM-2025 at a specified point/time.
 */
int main(void)
{
    //////////////////////////////////////////////////////////////////////
    const char *cof = "WMM2025.COF"; /* NOAA/NCEI COF file */

    // WMM initialization (reset state, load COF, internal init)
    wmm_init(cof);

    double epoch;
    char name[32];
    char date[16];
    wmm_version(&epoch, name, sizeof(name), date, sizeof(date));

    printf("# --------------------------------------------------\n");
    printf("# The World Magnetic Model (WMM) for %7.2lf\n", epoch);
    printf("# Valid year range: (%-7.2lf - %-7.2lf):\n", epoch, epoch + 5.0);
    printf("# Name : %s\n", name);
    printf("# Release date : %s\n", date);
    printf("# --------------------------------------------------\n");

    //////////////////////////////////////////////////////////////////////
    // input
    double dlat = 80.0;  /* deg */
    double dlon = 0.0;  /* deg */
    double alt  = 0.0;  /* km above WGS84 ellipsoid */
    double decimal_year = 0.0;

    //2025/01/01(YYYY/MM/DD) ==> 2025.0
    dec_year(2025, 1, 1, 0, 0, &decimal_year);

    printf(">>   Latitude(deg.) : %lf\n", dlat);
    printf(">>  Longitude(deg.) : %lf\n", dlon);
    printf(">>     Altitude(km) : %lf (WGS84 Ellipsoid)\n", alt);
    printf(">>     Decimal Yaer : %lf\n", decimal_year);

    //////////////////////////////////////////////////////////////////////
    // calculate with changes per year
    double decl, incl, ti1, gv;
    double x1, y1, z1, h1;
    double ax, ay, az, ah, af, ainc, adec;

    wmm_point_change(alt, dlat, dlon, decimal_year,
         &decl, &incl, &ti1, &gv,
         &x1, &y1, &z1, &h1,
         &ax, &ay, &az, &ah, &af,
         &adec, &ainc);

    MAGtype_GeoMagneticElements mag;
    memset(&mag, 0, sizeof(mag));
    mag.Decl = decl; mag.Incl = incl; mag.F = ti1; mag.GV = gv;
    mag.X = x1; mag.Y = y1; mag.Z = z1; mag.H = h1;

    mag.Decldot = adec; mag.Incldot = ainc; mag.Fdot = af;
    mag.Hdot = ah;     mag.Xdot = ax;       mag.Ydot = ay; mag.Zdot = az;
    mag.GVdot = mag.Decldot;

    //////////////////////////////////////////////////////////////////////
    // Uncertainty
    double ux, uy, uz, uh, uf, uinc, udec;
    wmm_error(h1, &uf, &uh, &ux, &uy, &uz, &uinc, &udec);

    MAGtype_GeoMagneticElements mag_error;
    memset(&mag_error, 0, sizeof(mag_error));
    mag_error.Decldot = udec;
    mag_error.Incldot = uinc;
    mag_error.F = uf; mag_error.H = uh; mag_error.X = ux; mag_error.Y = uy; mag_error.Z = uz;

    //////////////////////////////////////////////////////////////////////
    // Warning flags
    int warn_P = 0;
    int warn_H = 0;           double warn_H_val = 99999.0;
    int warn_H_strong = 0;    double warn_H_strong_val = 99999.0;

    if (h1 < 1000.0) { warn_H_strong = 1; warn_H_strong_val = h1; }
    else if (h1 < 5000.0) { warn_H = 1; warn_H_val = h1; }

    if (h1 < 100.0) { /* magnetic pole */
        decl = NAN;
        mag.Decl = decl;
    }

    //////////////////////////////////////////////////////////////////////
    wmm_print_results(dlat, dlon, alt, decimal_year,
                      &mag, &mag_error,
                      warn_H, warn_H_val,
                      warn_H_strong, warn_H_strong_val,
                      warn_P);

    // WMM initialization (reset state, load COF, internal init)
    wmm_init(cof);

    /* Run unit tests: tolerance ±0.1 */
    printf("# --------------------------------------------------\n");
    const char *tv = "WMM2025_TEST_VALUES.txt";
    run_wmm_tests(tv, 0.1);
}
