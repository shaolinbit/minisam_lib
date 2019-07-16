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
 * Find the full path to an example dataset。  The name
 * may be specified with or without a file extension - if no extension is
 * given, this function first looks for the .graph extension, then .txt.
 * @return The full path and filename to the requested dataset.
 * @throw std::invalid_argument if no matching file could be found using the
 * search process described above.
 */
std::string findExampleDataFile(const std::string& name);



#endif // DATASET_H
