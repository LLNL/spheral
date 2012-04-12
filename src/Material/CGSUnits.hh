//---------------------------------Spheral++----------------------------------//
// CGSUnits -- The base for the CGS unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_CGSUnits__
#define __Spheral_CGSUnits__

#include "PhysicalConstants.hh"

namespace Spheral {
namespace Material {

class CGSUnits: public PhysicalConstants<CGSUnits> {
public:
  //--------------------------- Public Interface ---------------------------//
  CGSUnits();
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
    class CGSUnits;
  }
}

#endif

  
