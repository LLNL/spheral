//---------------------------------Spheral++----------------------------------//
// SolarUnits -- The base for the Solar unit system.
//
// Created by JMO, Tue Mar 20 10:35:56 PDT 2012
//----------------------------------------------------------------------------//

#ifndef SolarUnits_HH
#define SolarUnits_HH

namespace Spheral {
namespace Material {

class SolarUnits {

public:
  //--------------------------- Public Interface ---------------------------//
  static const double unitLm;
  static const double unitMkg;
  static const double unitTsec;
};
}
}

#else
// Forward declaration.
namespace Spheral {
  namespace Material {
    class SolarUnits;
  }
}

#endif

  
