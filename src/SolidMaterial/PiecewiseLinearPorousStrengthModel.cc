//---------------------------------Spheral++----------------------------------//
// PiecewiseLinearPorousStrengthModel
//
// Shear modulus and yield are treated as piecewise linear functions of the
// the distension.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "PiecewiseLinearPorousStrengthModel.hh"
#include "StrengthModel.hh"
#include "Field/Field.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PiecewiseLinearPorousStrengthModel<Dimension>::
PiecewiseLinearPorousStrengthModel(const StrengthModel<Dimension>& solidStrength,
                                   const std::vector<Scalar> porosityAbscissa,
                                   const std::vector<Scalar> shearModulusRatios,
                                   const std::vector<Scalar> yieldStrengthRatios):
  PorousStrengthModel<Dimension>(solidStrength),
  mPorosityAbscissa(porosityAbscissa),
  mShearModulusRatios(shearModulusRatios),
  mYieldStrengthRatios(yieldStrengthRatios) {

    // make sure we're interpolating the full range
    if(mPorosityAbscissa.back() < 1.0){
        mPorosityAbscissa.push_back(1.0);
        mShearModulusRatios.push_back(0.0);
        mYieldStrengthRatios.push_back(0.0);
    }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PiecewiseLinearPorousStrengthModel<Dimension>::
~PiecewiseLinearPorousStrengthModel() {
}

//------------------------------------------------------------------------------
// Set the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PiecewiseLinearPorousStrengthModel<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& density,
             const Field<Dimension, Scalar>& specificThermalEnergy,
             const Field<Dimension, Scalar>& pressure,
             const Field<Dimension, SymTensor>& damage) const {
  REQUIRE(density.nodeListPtr() == shearModulus.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == shearModulus.nodeListPtr());
  REQUIRE(pressure.nodeListPtr() == shearModulus.nodeListPtr());
  REQUIRE(damage.nodeListPtr() == shearModulus.nodeListPtr());

  // The base model sets the solid (compacted) value.
  const Field<Dimension, Scalar>& alpha = this->alpha();
  const Field<Dimension, Scalar>  rhoS = alpha*density;
  const StrengthModel<Dimension>& strengthModel = this->solidStrength();
  
  strengthModel.shearModulus(shearModulus, rhoS, specificThermalEnergy, pressure, damage);

  // Now apply the porosity modifier.
  const unsigned n = shearModulus.numInternalElements();
  for (unsigned i = 0; i != n; ++i) {
    shearModulus(i) *= this->getPorousShearModulusRatio(alpha[i]);
  }
}

//------------------------------------------------------------------------------
// Set the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PiecewiseLinearPorousStrengthModel<Dimension>::
yieldStrength(Field<Dimension, Scalar>& yieldStrength,
              const Field<Dimension, Scalar>& density,
              const Field<Dimension, Scalar>& specificThermalEnergy,
              const Field<Dimension, Scalar>& pressure,
              const Field<Dimension, Scalar>& plasticStrain,
              const Field<Dimension, Scalar>& plasticStrainRate,
              const Field<Dimension, SymTensor>& damage) const {

  REQUIRE(density.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(pressure.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(damage.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(plasticStrain.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(plasticStrainRate.nodeListPtr() == yieldStrength.nodeListPtr());

  // The base model sets the solid (compacted) value.
  const Field<Dimension, Scalar>& alpha = this->alpha();
  const Field<Dimension, Scalar>  rhoS = alpha*density;
  const StrengthModel<Dimension>& strengthModel = this->solidStrength();
  
  strengthModel.yieldStrength(yieldStrength, rhoS, specificThermalEnergy, pressure, plasticStrain, plasticStrainRate, damage);

  // Now apply the porosity modifier.
  const unsigned n = yieldStrength.numInternalElements();
  for (unsigned i = 0; i != n; ++i) {
    yieldStrength(i) *= this->getPorousYieldStrengthRatio(alpha[i]);
  }
}

//------------------------------------------------------------------------------
// linear interpolation in phi
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PiecewiseLinearPorousStrengthModel<Dimension>::
getPorousYieldStrengthRatio(const Scalar alphai) const {

  // get our bin index  
  const unsigned numBins = mPorosityAbscissa.size()-1;

  Scalar YRatioInterp = 1.0;

  const Scalar phii = 1.0 - 1.0/alphai;
  
  for (unsigned i = 0; i != numBins; ++i) {
    const Scalar phi0 = mPorosityAbscissa[i];
    const Scalar phi1 = mPorosityAbscissa[i+1];
    if (phii >= phi0  and phii < phi1){
        const Scalar x = (phii - phi0)*safeInv(phi1-phi0);
        YRatioInterp = x*(mYieldStrengthRatios[i+1]-mYieldStrengthRatios[i])+mYieldStrengthRatios[i];
    }    
  }
  return YRatioInterp;
}

//------------------------------------------------------------------------------
// linear interpolation in phi
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PiecewiseLinearPorousStrengthModel<Dimension>::
getPorousShearModulusRatio(const Scalar alphai) const {

  // get our bin index  
  const unsigned numBins = mPorosityAbscissa.size()-1;

  Scalar muRatioInterp = 1.0;

  const Scalar phii = 1.0 - 1.0/alphai;
  
  for (unsigned i = 0; i != numBins; ++i) {
    const Scalar phi0 = mPorosityAbscissa[i];
    const Scalar phi1 = mPorosityAbscissa[i+1];
    if (phii >= phi0  and phii < phi1){
        const Scalar x = (phii - phi0)*safeInv(phi1-phi0);
        muRatioInterp = x*(mShearModulusRatios[i+1]-mShearModulusRatios[i])+mShearModulusRatios[i];
    }    
  }
  return muRatioInterp;
}

}
