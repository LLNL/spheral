//---------------------------------Spheral++----------------------------------//
// ConstantStrength -- An implentation of StrengthModel returning constant
// values for the shear modulus and yield strength.
//
// Created by JMO, Mon Jul 24 10:34:40 PDT 2006
//----------------------------------------------------------------------------//
#ifndef __Spheral_ConstantStrength_hh__
#define __Spheral_ConstantStrength_hh__

#include "StrengthModel.hh"
#include "SolidEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {

template<typename Dimension>
class ConstantStrength: public StrengthModel<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;

  // Constructors, destructor.
  ConstantStrength(const double mu0,
                   const double Y0);
  ConstantStrength(const double mu0,
                   const double Y0,
                   const SolidEquationOfState<Dimension>& eos);
  virtual ~ConstantStrength() {};

  // The generic interface we require all strength models to provide.
  virtual void shearModulus(Field<Dimension, Scalar>& shearModulus,
                            const Field<Dimension, Scalar>& density,
                            const Field<Dimension, Scalar>& specificThermalEnergy,
                            const Field<Dimension, Scalar>& pressure) const;

  virtual void yieldStrength(Field<Dimension, Scalar>& yieldStrength,
                             const Field<Dimension, Scalar>& density,
                             const Field<Dimension, Scalar>& specificThermalEnergy,
                             const Field<Dimension, Scalar>& pressure,
                             const Field<Dimension, Scalar>& plasticStrain,
                             const Field<Dimension, Scalar>& plasticStrainRate) const;

  // Read only access to the parameters.
  double mu0() const;
  double Y0() const; 

private:
  //--------------------------- Private Interface ---------------------------//
  // The values for shear modulus and yield strength.
  double mShearModulus0;
  double mYieldStrength0;
  const SolidEquationOfState<Dimension>* mEOSptr;

  // No copying or assignment.
  ConstantStrength(const ConstantStrength&);
  ConstantStrength& operator=(const ConstantStrength&);
};

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ConstantStrength<Dimension>::
ConstantStrength(const double mu0,
                 const double Y0):
  StrengthModel<Dimension>(),
  mShearModulus0(mu0),
  mYieldStrength0(Y0),
  mEOSptr(0) {
}

//------------------------------------------------------------------------------
// Constructor with SolidEquationOfState.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
ConstantStrength<Dimension>::
ConstantStrength(const double mu0,
                 const double Y0,
                 const SolidEquationOfState<Dimension>& eos):
  StrengthModel<Dimension>(),
  mShearModulus0(mu0),
  mYieldStrength0(Y0),
  mEOSptr(&eos) {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
ConstantStrength<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& /*density*/,
             const Field<Dimension, Scalar>& /*specificThermalEnergy*/,
             const Field<Dimension, Scalar>& /*pressure*/) const {
  shearModulus = mShearModulus0;
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
ConstantStrength<Dimension>::
yieldStrength(Field<Dimension, Scalar>& yieldStrength,
              const Field<Dimension, Scalar>& density,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/,
              const Field<Dimension, Scalar>& /*pressure*/,
              const Field<Dimension, Scalar>& /*plasticStrain*/,
              const Field<Dimension, Scalar>& /*plasticStrainRate*/) const {
  if (mEOSptr != 0 and
      density/(mEOSptr->referenceDensity()) < mEOSptr->etamin()) {
    yieldStrength = 0.0;
  } else {
    yieldStrength = mYieldStrength0;
  }
}

//------------------------------------------------------------------------------
// Return the input shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
ConstantStrength<Dimension>::
mu0() const {
  return mShearModulus0;
}

//------------------------------------------------------------------------------
// Return the input yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
ConstantStrength<Dimension>::
Y0() const {
  return mYieldStrength0;
}

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ConstantStrength;
}

#endif

