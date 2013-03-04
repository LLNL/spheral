//---------------------------------Spheral++----------------------------------//
// PorousStrengthModel
//
// An implementation of strain-alpha porosity model described in
// Wunnemann, Collins, & Melosh, Icarus, 180, 514-527 (2006)
// "A strain-based porosity model for use in hydrocode simulations of impacts
//  and implications for transient crater growth in porous targets"
//
// This model assumes you will provide a solid EOS which will be modified.
// The underlying actualy solid EOS should provide the reference density, which
// will be treated here as the compacted true solid reference density.
//
// Note this model introduces a new state variable, the distention (alpha), which
// the pressure now depends on.  This implies our usual definition of P(rho, eps)
// now becomes P(rho, eps, alpha).  Our EOS interface does not recognize this
// this parameter, so we store alpha locally and only allow Field updates of the
// pressure (forbidding the single value P lookup the EOS usually allows).
//
// Created by JMO, Thu Sep 13 15:50:28 PDT 2012
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorousStrengthModel_hh__
#define __Spheral_PorousStrengthModel_hh__

#include "StrengthModel.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class PorousStrengthModel: public StrengthModel<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;

  // Constructors, destructor.
  PorousStrengthModel(const StrengthModel<Dimension>& solidStrength);
  virtual ~PorousStrengthModel();

  // We require the Field only interface!
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

  // Forbid the non-Field calls.
  virtual double shearModulus(const double density,
                              const double specificThermalEnergy,
                              const double pressure) const { VERIFY2("PorousStrengthModel forbids non-Field shearModulus", false); return 0.0; }

  virtual double yieldStrength(const double density,
                               const double specificThermalEnergy,
                               const double pressure,
                               const double plasticStrain,
                               const double plasticStrainRate) const { VERIFY2("PorousStrengthModel forbids non-Field yieldStrength", false); return 0.0; }

  // Access the material parameters.
  const StrengthModel<Dimension>& solidStrength() const;
  const FieldSpace::Field<Dimension, Scalar>& alpha() const;
  void alpha(const FieldSpace::Field<Dimension, Scalar>& x);

private:
  //--------------------------- Private Interface ---------------------------//
  const StrengthModel<Dimension>& mSolidStrength;
  const FieldSpace::Field<Dimension, Scalar>* mAlphaPtr;

  // No copying or assignment.
  PorousStrengthModel(const PorousStrengthModel&);
  PorousStrengthModel& operator=(const PorousStrengthModel&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class PorousStrengthModel;
  }
}

#endif

