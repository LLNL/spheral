//---------------------------------Spheral++----------------------------------//
// MKSUnits -- The base for the MKS unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_MKSUnits__
#define __Spheral_MKSUnits__

#include "PhysicalConstants.hh"

namespace Spheral {
namespace Material {

class MKSUnits: public PhysicalConstants<MKSUnits> {

public:
  //--------------------------- Public Interface ---------------------------//
  MKSUnits();
  static const double unitLm;
  static const double unitMkg;
  static const double unitTsec;
};
}
}

#else
// Forward declaration.
namespace Spheral {
  namespace Kernel {
    class MKSUnits;
  }
}

#endif

  
