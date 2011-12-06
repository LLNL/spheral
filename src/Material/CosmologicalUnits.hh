//---------------------------------Spheral++----------------------------------//
// CosmologicalUnits -- The base for the Cosmological unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//

#ifndef CosmologicalUnits_HH
#define CosmologicalUnits_HH

namespace Spheral {
namespace Material {

class CosmologicalUnits {

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
    class CosmologicalUnits;
  }
}

#endif

  
