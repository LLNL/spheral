namespace Spheral {
//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MFVHydroBase<Dimension>::
DvolumeDt() const {
  return mDvolumeDt;
}
}