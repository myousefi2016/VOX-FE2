/*
 * Vox-FE (Version 2): Voxel-based Finite Element Modelling
 *
 * Interface developed as a ParaView plugin by
 * Richard Holbrey and Michael J. Fagan
 *
 * Dept of Mechanical Engineering
 * University of Hull
 * Cottingham Road
 * HU6 7RX
 *
 * Vox-FE (Version 1) was created by Andreas Bitternas, George Sisias, Jia Liu and others.
 */

#ifndef VOXFE_CONNECTED_COMPONENTS_IMAGE_DEFINES_H
#define VOXFE_CONNECTED_COMPONENTS_IMAGE_DEFINES_H


#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <limits>
#include <vector>
using std::cout;
using std::cerr;
using std::numeric_limits;
using std::string;
using std::set;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::fstream;
using std::endl;

#include <stdlib.h>

#if defined(_WIN32) || defined(_WIN64)
  #include <direct.h>
  #ifndef chdir
  #define chdir _chdir
  #endif
#else
  #include <unistd.h>
#endif

#define voxfeCCImageAppendString  "_v.mhd"
#define voxfeCCModelAppendString  "_v.model"
#define voxfeCCScriptAppendString "_v.script"
#define voxfeCCLogFileAppendString "_v.log"

void voxfeSplitFilename (const string& str, string& path, string& file, const char* sep )  {

  //std::cout << "Splitting: " << str << '\n';
  unsigned found = str.find_last_of(sep);

  //std::cout << " path: " << str.substr(0,found) << '\n';
  path = str.substr(0,found);

  //std::cout << " file: " << str.substr(found+1) << '\n';
  file = str.substr(found+1);
}


//a sstring comparator
struct voxfeStringLessThan {
    bool operator() (const string& lhs, const string& rhs) const {
      int result = lhs.compare(rhs);
      if( result < 0 ) return true;
      else if( result > 0 ) return false;

      return 0; //if equal
    }
};



#endif

