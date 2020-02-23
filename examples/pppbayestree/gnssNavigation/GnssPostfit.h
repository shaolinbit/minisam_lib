/**
 *  @file   GNSSPostfit.h
 *  @author Ryan
 *  @brief  Header file for GNSS postfit analysis
 **/

#pragma once


#include "../gnssNavigation/GnssData.h"
#include "../gnssNavigation/GnssTools.h"


#ifndef FOREACH_HPP
#define FOREACH_HPP
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#endif

#include "boost/foreach.hpp"

using namespace std;
using namespace boost;

#define foreach BOOST_FOREACH



/// Function to calculate gnss postfit residuals
vector<double> getResiduals( minivector &nomXYZ, std::map<int,minivector> &results,const vector<rnxData>& data);

/// write residuals to text file
void writeResiduals(const vector<double>& postfitResiduals, string outputFile,const vector<string>& switchIndex );

// iterate over residual vector to mark outliers.
vector<int> markResiduals(const vector<double>& postfitResdiuals, double threshold );


