//---------------------------------Spheral++----------------------------------//
// MKSUnits -- The base for the MKS unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//

#ifndef MKSUnits_HH
#define MKSUnits_HH

namespace Spheral {
namespace Material {

class MKSUnits {

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
  namespace Kernel {
    class MKSUnits;
  }
}

#endif

  
