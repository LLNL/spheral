//---------------------------------Spheral++----------------------------------//
// HydroFieldNames -- A collection of standard Field names for the hydro 
// physics package.
//----------------------------------------------------------------------------//
#ifndef _Spheral_GSPHFieldNames_
#define _Spheral_GSPHFieldNames_

#include <string>

namespace Spheral {

struct GSPHFieldNames {
  static const std::string mass;
};

}

#else

namespace Spheral {
  struct GSPHFieldNames;
}

#endif
