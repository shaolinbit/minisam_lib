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
#include "minisam/base/Matrix.h"
#include "minisam/base/MatCal.h"

using namespace std;

//// rotate from ECI to ECEF
Eigen::Vector3d inertialToECEF( const Eigen::Vector3d& inertialPosition, const double t, const double t0);

//// Generate rotation matrix from Earth-to-Navigation frame
Eigen::MatrixXd earthToNavTrans( const Eigen::Vector3d& ECEFxyz);

//// Compute mapping from meas. to states
Eigen::VectorXd obsMap(const Eigen::Vector3d& p1, const Eigen::Vector3d& q, const int& Trop = 1);

//// Compute mapping from meas. to states
Eigen::VectorXd  obsMapNED(const Eigen::Vector3d& p1, const Eigen::Vector3d& q, const int& Trop = 1);


//// Extract PRN Vector from GNSS data structure
Eigen::VectorXi getPRN(const Eigen::MatrixXd& p);

//// See if current PRN value was present at previous epoch
bool checkPRN(const Eigen::VectorXi& p, const int& n);

/// computer delta pseudorange observables
double deltaObs(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const double& pseudorange);

/// compute the delta troposphere correction
double deltaTrop(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);

//// Convert from WGS-84 ECEF coordinated to local-level-tangent (ENU) coordinates
////
//// REF :: Groves, Paul. Principles of GNSS, Inertial, and Multisensor Integrated
////        Navigation Systems. Artech House, 2008
Eigen::Vector3d xyz2enu(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);

//// Convert WGS-84 ECEF coordinates to LLH
////
//// REF :: Groves, Paul. Principles of GNSS, Inertial, and Multisensor Integrated
////        Navigation Systems. Artech House, 2008
Eigen::Vector3d xyz2llh(const Eigen::Vector3d& p1);

//// Convert ENU coordinates to ECEF
////
//// REF :: Groves, Paul. Principles of GNSS, Inertial, and Multisensor Integrated
////        Navigation Systems. Artech House, 2008
Eigen::Vector3d enu2xyz(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);

/// Convert NED Local Frame to ENU Local Frame
Eigen::Vector3d ned2enu(const Eigen::Vector3d& p1);

//// Computer Elevation of Satellite from Receiver
double calcEl(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2);

//// Compute elevation angle given a NED position vector.
double calcElNed(const Eigen::Vector3d& p1);

//// Computer elevation angle dependant weighting.
double elDepWeight(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, double measWeight);

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
double tropDry(const Eigen::Vector3d& p1);

//// time difference carrier-phase observations
double dopplerObs(const Eigen::Vector3d& p1, double tdcp1, const Eigen::Vector3d& p2, double tdcp2);

