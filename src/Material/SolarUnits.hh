//---------------------------------Spheral++----------------------------------//
// SolarUnits -- The base for the Solar unit system.
//
// Created by JMO, Tue Mar 20 10:35:56 PDT 2012
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolarUnits__
#define __Spheral_SolarUnits__

#include "PhysicalConstants.hh"

namespace Spheral {
namespace Material {

class SolarUnits: public PhysicalConstants<SolarUnits> {
public:
  //--------------------------- Public Interface ---------------------------//
  SolarUnits();
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

  
