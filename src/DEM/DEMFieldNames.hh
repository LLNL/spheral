//---------------------------------Spheral++----------------------------------//
// DEMFieldNames -- A collection of standard Field names for the DEM 
// physics package.
//
//----------------------------------------------------------------------------//
#ifndef _Spheral_DEMFieldNames_
#define _Spheral_DEMFieldNames_

#include <string>

namespace Spheral {

struct DEMFieldNames {
  static const std::string particleRadius;
  static const std::string angularVelocity;
  static const std::string uniqueIndices;
  static const std::string isActiveContact;
  static const std::string neighborIndices;
  static const std::string shearDisplacement;
  static const std::string equilibriumOverlap;
};

}

#else

namespace Spheral {
  struct DEMFieldNames;
}

#endif
