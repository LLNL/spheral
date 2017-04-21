namespace Spheral {
namespace CRKSPHSpace {

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

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SolidCRKSPHHydroBase<Dimension>::
Adamage() const {
  return mAdamage;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SolidCRKSPHHydroBase<Dimension>::
Bdamage() const {
  return mBdamage;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SolidCRKSPHHydroBase<Dimension>::
Cdamage() const {
  return mCdamage;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SolidCRKSPHHydroBase<Dimension>::
gradAdamage() const {
  return mGradAdamage;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SolidCRKSPHHydroBase<Dimension>::
gradBdamage() const {
  return mGradBdamage;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>&
SolidCRKSPHHydroBase<Dimension>::
gradCdamage() const {
  return mGradCdamage;
}

}
}
