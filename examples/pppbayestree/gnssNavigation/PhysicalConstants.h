/**
 * @file   PhysicalConstants.h
 * @brief  House values that are repeatedly used
 * @author Ryan Watson
 */

#pragma once

//namespace gtsam {

const int GNSS_speedOfLight = 299792458;   // [m/s]
const int GNSS_gravity = 9.80665;   // [m/s*s]
const double GNSS_L1 = 1575.42e+6;   // GPS L1 freq [HZ]
const double GNSS_L2 = 1227.6e+6;   // GPS L2 freq [Hz]
const double GNSS_L5 = 1176.45e+6;   // GPS L5 freq [Hz]
const double GNSS_semiMajor = 6378137.0;  // Earth's Semi-Major Axis [m]
const double GNSS_semiMinor = 6356752.3142;    // Earth's Semi-Minor Axis [m]
const double GNSS_earthRot = 7.292115e-5;     // Earth's rotation rate [rad/sec]
const int GNSS_stdPressure = 1013;   // Std Atmosphere Pressure [mbar]
const double GNSS_stdTemp = 288.15;    // Std Atmosphere Temp. [Kelvin]

//}
