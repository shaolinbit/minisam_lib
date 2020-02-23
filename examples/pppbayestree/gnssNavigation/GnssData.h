/**
 * @file   GnssData.h
 * @brief  Tools required to read/write GNSS data
 * @author Ryan Watson
 */


#pragma once

#include "../gnssNavigation/GnssTools.h"
#include "minisam/miniblas/minivector_double.h"
#include "minisam/miniblas/minimatrix_double.h"


#include "boost/foreach.hpp"
#include "boost/tuple/tuple.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
using namespace boost;

#define foreach BOOST_FOREACH



/// Read GNSS data in the rnxToGtsam.cpp format
/// Data = { Week, Sow, Epoch, SVN, SatXYZ, Rho, P.C., L.C., Break_Flag}
typedef boost::tuple<double, int, int, minivector, double, double, double, int> rnxData;
vector<rnxData> readGNSS(const std::string& fileLoc,const std::string& filepath=NULL);

/// Write GNSS states to text file
void writeStates(std::map<int,minivector> &results, string outputFile);

/// Write Pos. solution in ENU co-ordinate frame
void writeNavFrame(std::map<int,minivector>  &results, minivector &nom, string outputFile);

/// Write Pos. solution in ECEF co-ordinate frame
void writeEarthFrame(std::map<int,minivector>  &results, minivector &nom, string outputFile);

/// Write switch states to text file
void writeSwitches( std::map<int,minivector>  &results, string outputFile, vector<string> switchIndex);

/// Write switch states to text file
void writeAmbiguity( std::map<int,minivector>  &results, string outputFile, vector<string> satIndex);


