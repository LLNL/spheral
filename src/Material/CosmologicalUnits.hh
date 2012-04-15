//---------------------------------Spheral++----------------------------------//
// CosmologicalUnits -- The base for the Cosmological unit system.
//
// Created by JMO, Fri Mar 31 17:07:41 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_CosmologicalUnits__
#define __Spheral_CosmologicalUnits__

#include "PhysicalConstants.hh"

namespace Spheral {
namespace Material {

class CosmologicalUnits: public PhysicalConstants<CosmologicalUnits> {
public:
  //--------------------------- Public Interface ---------------------------//
  CosmologicalUnits();
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

  
