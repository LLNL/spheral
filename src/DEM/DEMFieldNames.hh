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
};

}

#else

namespace Spheral {
  struct DEMFieldNames;
}

#endif
