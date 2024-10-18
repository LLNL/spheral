namespace Spheral {

//------------------------------------------------------------------------------
// Control whether allow damaged material to have stress relieved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidSPHHydroBase<Dimension>::
damageRelieveRubble() const {
  return mDamageRelieveRubble;
}

template<typename Dimension>
inline
void
SolidSPHHydroBase<Dimension>::
damageRelieveRubble(bool x) {
  mDamageRelieveRubble = x;
}

//------------------------------------------------------------------------------
// Do we allow damaged material to have strength?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidSPHHydroBase<Dimension>::
strengthInDamage() const {
  return mStrengthInDamage;
}

template<typename Dimension>
inline
void
SolidSPHHydroBase<Dimension>::
strengthInDamage(bool x) {
  mStrengthInDamage = x;
}

//------------------------------------------------------------------------------
// Kernel for velocity gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
SolidSPHHydroBase<Dimension>::
GradKernel() const {
  return mGradKernel;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SolidSPHHydroBase<Dimension>::
DdeviatoricStressDt() const {
  return mDdeviatoricStressDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidSPHHydroBase<Dimension>::
bulkModulus() const {
  return mBulkModulus;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidSPHHydroBase<Dimension>::
shearModulus() const {
  return mShearModulus;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidSPHHydroBase<Dimension>::
yieldStrength() const {
  return mYieldStrength;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidSPHHydroBase<Dimension>::
plasticStrain0() const {
  return mPlasticStrain0;
}

}
