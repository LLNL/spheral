//---------------------------------Spheral++----------------------------------//
// CGSUnits -- The base for the CGS unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//

#ifndef CGSUnits_HH
#define CGSUnits_HH

namespace Spheral {
namespace Material {

class CGSUnits {

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
    class CGSUnits;
  }
}

#endif

  
