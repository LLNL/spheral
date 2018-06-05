namespace Spheral {
namespace CRKSPHSpace {

//------------------------------------------------------------------------------
// Control whether allow damaged material to have stress relieved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SolidCRKSPHHydroBase<Dimension>::
damageRelieveRubble() const {
  return mDamageRelieveRubble;
}

template<typename Dimension>
inline
void
SolidCRKSPHHydroBase<Dimension>::
damageRelieveRubble(const bool x) {
  mDamageRelieveRubble = x;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SolidCRKSPHHydroBase<Dimension>::
DdeviatoricStressDt() const {
  return mDdeviatoricStressDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SolidCRKSPHHydroBase<Dimension>::
bulkModulus() const {
  return mBulkModulus;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SolidCRKSPHHydroBase<Dimension>::
shearModulus() const {
  return mShearModulus;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SolidCRKSPHHydroBase<Dimension>::
yieldStrength() const {
  return mYieldStrength;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SolidCRKSPHHydroBase<Dimension>::
plasticStrain0() const {
  return mPlasticStrain0;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SolidCRKSPHHydroBase<Dimension>::
Hfield0() const {
  return mHfield0;
}

}
}
