namespace Spheral {
//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
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
}