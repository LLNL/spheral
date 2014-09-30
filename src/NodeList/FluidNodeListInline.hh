#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "SmoothingScaleBase.hh"

namespace Spheral {
namespace NodeSpace {

//------------------------------------------------------------------------------
// Mass density per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldSpace::Field<Dimension, typename Dimension::Scalar>&
FluidNodeList<Dimension>::massDensity() {
  REQUIRE(mMassDensity.nodeListPtr() == this);
  return mMassDensity;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
FluidNodeList<Dimension>::massDensity() const {
  REQUIRE(mMassDensity.nodeListPtr() == this);
  return mMassDensity;
}

//------------------------------------------------------------------------------
// Specific thermal energy per node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldSpace::Field<Dimension, typename Dimension::Scalar>&
FluidNodeList<Dimension>::specificThermalEnergy() {
  REQUIRE(mSpecificThermalEnergy.nodeListPtr() == this);
  return mSpecificThermalEnergy;
}

template<typename Dimension>
inline
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
FluidNodeList<Dimension>::specificThermalEnergy() const {
  REQUIRE(mSpecificThermalEnergy.nodeListPtr() == this);
  return mSpecificThermalEnergy;
}

//------------------------------------------------------------------------------
// Access the equation of state.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const Material::EquationOfState<Dimension>&
FluidNodeList<Dimension>::equationOfState() const {
  return *mEosPtr;
}

template<typename Dimension>
inline
void
FluidNodeList<Dimension>::equationOfState(const Material::EquationOfState<Dimension>& eos) {
  mEosPtr = &eos;
}

//------------------------------------------------------------------------------
// Min/max mass densities when we're time integrating.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
FluidNodeList<Dimension>::
rhoMin() const {
  return mRhoMin;
}

template<typename Dimension>
inline
void
FluidNodeList<Dimension>::
rhoMin(const typename Dimension::Scalar x) {
  mRhoMin = x;
}

template<typename Dimension>
inline
typename Dimension::Scalar
FluidNodeList<Dimension>::
rhoMax() const {
  return mRhoMax;
}

template<typename Dimension>
inline
void
FluidNodeList<Dimension>::
rhoMax(const typename Dimension::Scalar x) {
  mRhoMax = x;
}

}
}
