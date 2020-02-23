#ifndef DATASET_H
#define DATASET_H


/**
 * @file dataset.h
 * @date
 * @author
 * @brief utility functions for loading datasets
 */

#pragma once
#include <string>
#include <utility> // for pair
#include <vector>
#include <iosfwd>


/**
 * Find the full path to an example dataset distributed with gtsam.  The name
 * may be specified with or without a file extension - if no extension is
 * given, this function first looks for the .graph extension, then .txt.  We
 * first check the gtsam source tree for the file, followed by the installed
 * example dataset location.  Both the source tree and installed locations
 * are obtained from CMake during compilation.
 * @return The full path and filename to the requested dataset.
 * @throw std::invalid_argument if no matching file could be found using the
 * search process described above.
 */
std::string findExampleDataFile(const std::string& name,const std::string& path=NULL);



#endif // DATASET_H
