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

#include "SolidEquationOfState.hh"
#include "DataOutput/registerWithRestart.hh"

#include <limits>

namespace Spheral {

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
  PorousEquationOfState(const EquationOfState<Dimension>& EOS);
  virtual ~PorousEquationOfState();

  //............................................................................
  // EOS methods.
  // We require any equation of state to define the following methods for Fields.
  virtual void setPressure(Field<Dimension, Scalar>& Pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const;

  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const;

  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setGammaField(Field<Dimension, Scalar>& gamma,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const;

  // Check if the underlying SolidEquationOfState is valid.
  virtual bool valid() const;
  //............................................................................

  // Access the material parameters.
  const EquationOfState<Dimension>& EOS() const;
  const SolidEquationOfState<Dimension>& solidEOS() const;  // Throws if the underlying EOS is not a SolidEquationOfState

  const Field<Dimension, Scalar>& alpha() const;
  void alpha(const Field<Dimension, Scalar>& x);

  const Field<Dimension, Scalar>& alpha0() const;
  void alpha0(const Field<Dimension, Scalar>& x);

  const Field<Dimension, Scalar>& c0() const;
  void c0(const Field<Dimension, Scalar>& x);

private:
  //--------------------------- Private Interface ---------------------------//
  const EquationOfState<Dimension>& mEOS;
  const Field<Dimension, Scalar> *mAlphaPtr, *mAlpha0Ptr, *mC0Ptr;

  // Disallow default constructor
  PorousEquationOfState();
};

}

#endif
