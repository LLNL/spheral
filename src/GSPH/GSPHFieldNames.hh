//---------------------------------Spheral++----------------------------------//
// HydroFieldNames -- A collection of standard Field names for the hydro 
// physics package.
//----------------------------------------------------------------------------//
#ifndef _Spheral_GSPHFieldNames_
#define _Spheral_GSPHFieldNames_

#include <string>

namespace Spheral {

struct GSPHFieldNames {
  static const std::string densityGradient;
  static const std::string pressureGradient;
  static const std::string deviatoricStressTensorGradient;
  static const std::string RiemannPressureGradient;
  static const std::string RiemannVelocityGradient;
  static const std::string RiemannDeviatoricStressTensorGradient;
};

}

#else

namespace Spheral {
  struct GSPHFieldNames;
}

#endif
