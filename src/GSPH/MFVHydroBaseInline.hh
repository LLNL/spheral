namespace Spheral {


//------------------------------------------------------------------------------
// set/get for mesh motion coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
MFVHydroBase<Dimension>::nodeMotionCoefficient() const {
  return mNodeMotionCoefficient;
}

template<typename Dimension>
inline
void
MFVHydroBase<Dimension>::
nodeMotionCoefficient(typename Dimension::Scalar x) {
  mNodeMotionCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get mesh motion type
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NodeMotionType
MFVHydroBase<Dimension>::
nodeMotionType() const {
  return mNodeMotionType;
}

template<typename Dimension>
inline
void
MFVHydroBase<Dimension>::
nodeMotionType(NodeMotionType x) {
  mNodeMotionType=x;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
MFVHydroBase<Dimension>::
nodalVelocity() const {
  return mNodalVelocity;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MFVHydroBase<Dimension>::
DmassDt() const {
  return mDmassDt;
}
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MFVHydroBase<Dimension>::
DthermalEnergyDt() const {
  return mDthermalEnergyDt;
}
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
MFVHydroBase<Dimension>::
DmomentumDt() const {
  return mDmomentumDt;
}
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MFVHydroBase<Dimension>::
DvolumeDt() const {
  return mDvolumeDt;
}
template<typename Dimension>
inline
const std::vector<typename Dimension::Scalar>&
MFVHydroBase<Dimension>::
pairMassFlux() const {
  return mPairMassFlux;
}
// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::SymTensor>&
// MFVHydroBase<Dimension>::
// HStretchTensor() const {
//   return mHStretchTensor;
// }
}