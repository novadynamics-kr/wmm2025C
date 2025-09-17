/*======================================================================*
 *  File        : wmm_mag_struct.h
 *  Brief       : World Magnetic Model (WMM-2025) public API (C)
 *  Author      : Nova Dynamics (novadynamics.kr@gmail.com)
 *  License     : MIT
 *
 *  Notes
 *  - Units:
 *      * latitude/longitude: degrees (geodetic; East/North positive)
 *      * altitude          : kilometers above WGS84 ellipsoid
 *      * time              : decimal year (e.g., 2025.0)
 *      * field outputs     : nT; angles in degrees
 *======================================================================*/

#ifndef WMM_MAG_STRUCT_H
#define WMM_MAG_STRUCT_H

/**
 * @brief Structure for storing computed geomagnetic elements and secular variation.
 * @details
 * - Each field stores either the WMM calculation result or the annual change (difference after 1 year).
 * - Units: angles (deg), magnetic field (nT), secular variation (nT/yr or deg/yr)
 */
typedef struct {
    double Decl;    /**< @brief Declination D (deg). Angle between magnetic north and true north, east-positive */
    double Incl;    /**< @brief Inclination I (deg). Angle between horizontal plane and magnetic field vector, down-positive */
    double F;       /**< @brief Total field intensity F (nT) */
    double H;       /**< @brief Horizontal component H (nT) = sqrt(X^2 + Y^2) */
    double X;       /**< @brief North component X (nT) = F cos(D) cos(I) */
    double Y;       /**< @brief East component Y (nT) = F sin(D) cos(I) */
    double Z;       /**< @brief Down component Z (nT) = F sin(I) */
    double GV;      /**< @brief Grid variation GV (deg, meaningful only at high latitudes) */
    double Decldot; /**< @brief Annual change of declination dD/dt (deg/yr) */
    double Incldot; /**< @brief Annual change of inclination dI/dt (deg/yr) */
    double Fdot;    /**< @brief Annual change of total intensity dF/dt (nT/yr) */
    double Hdot;    /**< @brief Annual change of horizontal component dH/dt (nT/yr) */
    double Xdot;    /**< @brief Annual change of north component dX/dt (nT/yr) */
    double Ydot;    /**< @brief Annual change of east component dY/dt (nT/yr) */
    double Zdot;    /**< @brief Annual change of down component dZ/dt (nT/yr) */
    double GVdot;   /**< @brief Annual change of GV dGV/dt (deg/yr) */
} MAGtype_GeoMagneticElements;

/**
 * @brief WMM internal calculation state structure.
 * @details
 * - Contains all internal variables and caches required for WMM calculation.
 * - Stores all intermediate calculation values such as Gauss coefficients read from COF file, Schmidt normalization, Legendre functions, longitude Fourier terms, etc.
 */
typedef struct {
    int   maxord;           /**< Maximum degree/order (usually 12 for WMM) */

    double c[13][13];       /**< Gauss coefficients g_nm (nT), main geomagnetic model */
    double cd[13][13];      /**< Gauss coefficient secular variation dg_nm/dt (nT/yr) */
    double tc[13][13];      /**< Time-corrected coefficients g_nm(t) = g_nm + dt*dg_nm */
    double dp[13][13];      /**< Derivative of Legendre function ∂P̄_n^m/∂θ (theta-derivative) */
    double snorm[169];      /**< Schmidt semi-normalization coefficients (index n+m*13) */
    double sp[13];          /**< sin(mλ): longitude-dependent Fourier terms */
    double cp[13];          /**< cos(mλ): longitude-dependent Fourier terms */
    double fn[13];          /**< n+1: recursion helper */
    double fm[13];          /**< m: recursion helper */
    double pp[13];          /**< Auxiliary array for pole special handling */
    double k[13][13];       /**< Legendre recursion coefficients k[m][n] */

    double pi;              /**< Pi value (M_PI) */
    double dtr;             /**< Degree-to-radian conversion constant (pi/180) */

    double a;               /**< WGS84 semi-major axis (km) */
    double b;               /**< WGS84 semi-minor axis (km) */
    double re;              /**< WMM reference radius (6371.2 km) */

    double a2;              /**< a^2: semi-major axis squared */
    double b2;              /**< b^2: semi-minor axis squared */
    double c2;              /**< a^2 - b^2 */
    double a4;              /**< a^4: semi-major axis to the fourth power */
    double b4;              /**< b^4: semi-minor axis to the fourth power */
    double c4;              /**< a^4 - b^4 */

    double epoch;           /**< Model reference epoch (e.g., 2025.0) */

    // Last evaluation point cache (for performance optimization)
    double otime;           /**< Last evaluation time (decimal year) */
    double oalt;            /**< Last evaluation altitude (km) */
    double olat;            /**< Last evaluation latitude (deg) */
    double olon;            /**< Last evaluation longitude (deg) */

    // Geocentric coordinate transformation result cache
    double ct;              /**< cos(theta) = cos(geocentric colatitude) */
    double st;              /**< sin(theta) */
    double r;               /**< Geocentric radius (km) */
    double ca;              /**< cos(geodetic→geocentric tilt) */
    double sa;              /**< sin(geodetic→geocentric tilt) */
} WMMState;

/**
 * @brief WMM version and metadata structure.
 * @details
 * - Stores model name, reference epoch, release date, etc. read from COF file header.
 */
typedef struct {
    double epoch;           /**< @brief Model valid start year (decimal year, e.g., 2025.0) */
    char   model_name[32];  /**< @brief Model name (e.g., "WMM-2025") */
    char   release_date[16];/**< @brief Release date "MM/DD/YYYY" (optional) */
} WMMVersion;

#endif /* WMM_MAG_STRUCT_H */