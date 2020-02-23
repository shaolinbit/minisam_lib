/**
 * @file   GnssTools.h
 * @brief  Tools required to process GNSS data -- (i.e. ECEF to ENU transformation)
 * @author Ryan Watson & Jason Gross
 */

#pragma once

#include "../gnssNavigation/PhysicalConstants.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include "minisam/mat/Matrix.h"
#include "minisam/mat/MatCal.h"

using namespace std;

//// rotate from ECI to ECEF
minivector inertialToECEF( const minivector& inertialPosition, const double t, const double t0);

//// Generate rotation matrix from Earth-to-Navigation frame
minimatrix earthToNavTrans( const minivector& ECEFxyz);

//// Compute mapping from meas. to states
minivector obsMap(const minivector& p1, const minivector& q, const int& Trop = 1);

//// Compute mapping from meas. to states
minivector  obsMapNED(const minivector& p1, const minivector& q, const int& Trop = 1);


//// Extract PRN Vector from GNSS data structure
mini_int_vector getPRN(const minimatrix& p);

//// See if current PRN value was present at previous epoch
bool checkPRN(const mini_int_vector& p, const int& n);

/// computer delta pseudorange observables
double deltaObs(const minivector& p1, const minivector& p2, const double& pseudorange);

/// compute the delta troposphere correction
double deltaTrop(const minivector& p1, const minivector& p2);

//// Convert from WGS-84 ECEF coordinated to local-level-tangent (ENU) coordinates
////
//// REF :: Groves, Paul. Principles of GNSS, Inertial, and Multisensor Integrated
////        Navigation Systems. Artech House, 2008
minivector xyz2enu(const minivector& p1, const minivector& p2);

//// Convert WGS-84 ECEF coordinates to LLH
////
//// REF :: Groves, Paul. Principles of GNSS, Inertial, and Multisensor Integrated
////        Navigation Systems. Artech House, 2008
minivector xyz2llh(const minivector& p1);

//// Convert ENU coordinates to ECEF
////
//// REF :: Groves, Paul. Principles of GNSS, Inertial, and Multisensor Integrated
////        Navigation Systems. Artech House, 2008
minivector enu2xyz(const minivector& p1, const minivector& p2);

/// Convert NED Local Frame to ENU Local Frame
minivector ned2enu(const minivector& p1);

//// Computer Elevation of Satellite from Receiver
double calcEl(const minivector& p1, const minivector& p2);

//// Compute elevation angle given a NED position vector.
double calcElNed(const minivector& p1);

//// Computer elevation angle dependant weighting.
double elDepWeight(const minivector& p1, const minivector& p2, double measWeight);

//// Elevation angle only troposphere mapping
////
//// REF :: Black, H. and Eisner, A., 1984. Correcting satellite Doppler data for
////        tropospheric effects. Journal of Geophysical Research.
double tropMap(const double& El);

//// Troposphere Model --- dry component only. Uses the Sass. model.
////
//// REF :: Saastamoinen, J. 1972, Atmospheric correction for the troposphere and
////        stratosphere in radio ranging of satellites, in The Use of Artificial
////        Satellites for Geodesy
double tropDry(const minivector& p1);

//// time difference carrier-phase observations
double dopplerObs(const minivector& p1, double tdcp1, const minivector& p2, double tdcp2);

