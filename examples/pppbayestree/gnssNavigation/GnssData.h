/**
 * @file   GnssData.h
 * @brief  Tools required to read/write GNSS data
 * @author Ryan Watson
 */


#pragma once

#include "../gnssNavigation/GnssTools.h"

#include "boost/foreach.hpp"
#include "boost/tuple/tuple.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
using namespace boost;

#define foreach BOOST_FOREACH



/// Read GNSS data in the format
/// Data = { Week, Sow, Epoch, SVN, SatXYZ, Rho, P.C., L.C., Break_Flag}
typedef boost::tuple<double, int, int, Eigen::Vector3d, double, double, double, int> rnxData;
vector<rnxData> readGNSS(const std::string& fileLoc);

/// Write GNSS states to text file
void writeStates(std::map<int,Eigen::VectorXd> &results, string outputFile);

/// Write Pos. solution in ENU co-ordinate frame
void writeNavFrame(std::map<int,Eigen::VectorXd> &results, Eigen::Vector3d &nom, string outputFile);

/// Write Pos. solution in ECEF co-ordinate frame
void writeEarthFrame(std::map<int,Eigen::VectorXd> &results, Eigen::Vector3d &nom, string outputFile);

/// Write switch states to text file
void writeSwitches( std::map<int,Eigen::VectorXd> &results, string outputFile, vector<string> switchIndex);

/// Write switch states to text file
void writeAmbiguity( std::map<int,Eigen::VectorXd> &results, string outputFile, vector<string> satIndex);


