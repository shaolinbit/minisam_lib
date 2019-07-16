/**
 * @file dataset.cpp
 * @date
 * @author
 * @brief utility functions for loading datasets
 */


#include "boost/filesystem/operations.hpp"

#include "boost/filesystem/path.hpp"

//#define BOOST_NO_CXX11_SCOPED_ENUMS
#include "boost/filesystem.hpp"
//#undef BOOST_NO_CXX11_SCOPED_ENUMS


#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace boost;
namespace fs = boost::filesystem;


string findExampleDataFile(const string& name)
{
  // Search source tree and installed location
  std::vector<string> rootsToSearch;
  rootsToSearch.push_back("pppbayestree\\data");//windows setting.
  rootsToSearch.push_back("pppbayestree/data");//ubuntu setting.

  // Search for filename as given, and with .graph and .txt extensions
  std::vector<string> namesToSearch;
  namesToSearch.push_back(name);
  namesToSearch.push_back(name + ".graph");
  namesToSearch.push_back(name + ".txt");
  namesToSearch.push_back(name + ".out");

  // Find first name that exists
  for(const fs::path& root: rootsToSearch) {
    for(const fs::path& name: namesToSearch) {
    //  cout<< (root / name).string()<<endl;
      if (fs::is_regular_file(root / name))
        return (root / name).string();
    }
  }

  // If we did not return already, then we did not find the file
  throw invalid_argument(
      "findExampleDataFile could not find a matching file in ../pppbayestree/data");
}

