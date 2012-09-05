//---------------------------------Spheral++----------------------------------//
// PorousEquationOfState
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
// Created by JMO, Fri Jun  1 16:16:26 PDT 2012
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorousEquationOfState__
#define __Spheral_PorousEquationOfState__

#include <limits>
#include "SolidEquationOfState.hh"
#include "DataOutput/registerWithRestart.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class PorousEquationOfState: 
    public SolidEquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  PorousEquationOfState(const SolidEquationOfState<Dimension>& solidEOS); // Solid EOS we're going to modify
  virtual ~PorousEquationOfState();

  //............................................................................
  // EOS methods.
  // We require any equation of state to define the following methods for Fields.
  virtual void setPressure(FieldSpace::Field<Dimension, Scalar>& Pressure,
                           const FieldSpace::Field<Dimension, Scalar>& massDensity,
                           const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setTemperature(FieldSpace::Field<Dimension, Scalar>& temperature,
                              const FieldSpace::Field<Dimension, Scalar>& massDensity,
                              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setSpecificThermalEnergy(FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                                        const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                        const FieldSpace::Field<Dimension, Scalar>& temperature) const;

  virtual void setSpecificHeat(FieldSpace::Field<Dimension, Scalar>& specificHeat,
                               const FieldSpace::Field<Dimension, Scalar>& massDensity,
                               const FieldSpace::Field<Dimension, Scalar>& temperature) const;

  virtual void setSoundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
                             const FieldSpace::Field<Dimension, Scalar>& massDensity,
                             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setGammaField(FieldSpace::Field<Dimension, Scalar>& gamma,
			     const FieldSpace::Field<Dimension, Scalar>& massDensity,
			     const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setBulkModulus(FieldSpace::Field<Dimension, Scalar>& bulkModulus,
			     const FieldSpace::Field<Dimension, Scalar>& massDensity,
			     const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  // Some of the following methods are disabled because we need to know alpha as well.
  virtual Scalar pressure(const Scalar massDensity,
                          const Scalar specificThermalEnergy) const { VERIFY2(false, "PorousEquationOfState does not support individual state calls."); }

  virtual Scalar temperature(const Scalar massDensity,
                             const Scalar specificThermalEnergy) const { VERIFY2(false, "PorousEquationOfState does not support individual state calls."); }

  virtual Scalar specificThermalEnergy(const Scalar massDensity,
                                       const Scalar temperature) const { VERIFY2(false, "PorousEquationOfState does not support individual state calls."); }

  virtual Scalar specificHeat(const Scalar massDensity,
                              const Scalar temperature) const { VERIFY2(false, "PorousEquationOfState does not support individual state calls."); }

  virtual Scalar soundSpeed(const Scalar massDensity,
                            const Scalar specificThermalEnergy) const { VERIFY2(false, "PorousEquationOfState does not support individual state calls."); }

  // Get the effective gamma (ratio of specific heats) for this eos.
  virtual Scalar gamma(const Scalar massDensity,
		       const Scalar specificThermalEnergy) const { VERIFY2(false, "PorousEquationOfState does not support individual state calls."); }

  // Get the bulk modulus.
  virtual Scalar bulkModulus(const Scalar massDensity,
                             const Scalar specificThermalEnergy) const { VERIFY2(false, "PorousEquationOfState does not support individual state calls."); }

  // Check if the underlying SolidEquationOfState is valid.
  virtual bool valid() const;
  //............................................................................

  // Access the material parameters.
  const SolidEquationOfState<Dimension>& solidEOS() const;
  const FieldSpace::Field<Dimension, Scalar>& alpha() const;
  void alpha(const FieldSpace::Field<Dimension, Scalar>& x);

private:
  //--------------------------- Private Interface ---------------------------//
  const SolidEquationOfState<Dimension>& mSolidEOS;
  const FieldSpace::Field<Dimension, Scalar>* mAlphaPtr;

  // Disallow default constructor
  PorousEquationOfState();
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class PorousEquationOfState;
  }
}

#endif
