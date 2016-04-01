//---------------------------------Spheral++----------------------------------//
// PhysicsEvolvingMaterialLibrary -- A base class to combine the Physics,
//  EquationOfState, and StrengthModel for more sophisticated material
//  libraries.
//
// This class implements three distinct Spheral interfaces:
//   EquationOfState : provides the ordinary Spheral EOS interface
//   StrengthModel   : also answers the Strength questions about shear modulus 
//                     and yield strength
//   Physics         : Geodyn also needs to advance it's own internal state, so
//                     this class implements the Physics interface to support
//                     that.
//
// Created by JMO, Thu Mar 31 16:11:26 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_PhysicsEvolvingMaterialLibrary_hh__
#define __Spheral_PhysicsEvolvingMaterialLibrary_hh__

#include "boost/multi_array.hpp"

#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Physics/Physics.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class PhysicsEvolvingMaterialLibrary: 
    public Material::EquationOfState<Dimension>,
    public SolidMaterial::StrengthModel<Dimension>,
    public PhysicsSpace::Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors.
  // We should add arguments to the constructor to specify the PhysicsEvolvingMaterialLibrary model
  PhysicsEvolvingMaterialLibrary(const Material::PhysicalConstants& constants,
                                 const double minimumPressure,
                                 const double maximumPressure,
                                 const Material::MaterialPressureMinType minPressureType);
  virtual ~PhysicsEvolvingMaterialLibrary();

private:
  //--------------------------- Private Interface ---------------------------//
  // Disallow default constructor
  PhysicsEvolvingMaterialLibrary();
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class PhysicsEvolvingMaterialLibrary;
  }
}

#endif
