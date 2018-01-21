//---------------------------------Spheral++----------------------------------//
// NullStrength -- An implentation of StrengthModel which mimics having
// no strength.
//
// Created by JMO, Thu Sep 16 23:48:47 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_NullStrength_hh__
#define __Spheral_NullStrength_hh__

#include "StrengthModel.hh"
#include "Field/Field.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class NullStrength: public StrengthModel<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;

  // Constructors, destructor.
  NullStrength() {};
  virtual ~NullStrength() {};

  // The generic interface we require all strength models to provide.
  virtual void shearModulus(FieldSpace::Field<Dimension, Scalar>& shearModulus,
                            const FieldSpace::Field<Dimension, Scalar>& density,
                            const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                            const FieldSpace::Field<Dimension, Scalar>& pressure) const;

  virtual void yieldStrength(FieldSpace::Field<Dimension, Scalar>& yieldStrength,
                             const FieldSpace::Field<Dimension, Scalar>& density,
                             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                             const FieldSpace::Field<Dimension, Scalar>& pressure,
                             const FieldSpace::Field<Dimension, Scalar>& plasticStrain,
                             const FieldSpace::Field<Dimension, Scalar>& plasticStrainRate) const;

private:
  //--------------------------- Private Interface ---------------------------//

  // No copying or assignment.
  NullStrength(const NullStrength&);
  NullStrength& operator=(const NullStrength&);
};

#ifndef __GCCXML__

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NullStrength<Dimension>::
shearModulus(FieldSpace::Field<Dimension, Scalar>& shearModulus,
             const FieldSpace::Field<Dimension, Scalar>& density,
             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
             const FieldSpace::Field<Dimension, Scalar>& pressure) const {
  shearModulus = 0.0;
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
NullStrength<Dimension>::
yieldStrength(FieldSpace::Field<Dimension, Scalar>& yieldStrength,
              const FieldSpace::Field<Dimension, Scalar>& density,
              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
              const FieldSpace::Field<Dimension, Scalar>& pressure,
              const FieldSpace::Field<Dimension, Scalar>& plasticStrain,
              const FieldSpace::Field<Dimension, Scalar>& plasticStrainRate) const {
  yieldStrength = 0.0;
}

#endif

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class NullStrength;
  }
}

#endif

