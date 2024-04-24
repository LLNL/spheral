namespace Spheral {

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
SPHSmoothingScale<Dimension>::
WT() const {
  return mWT;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHSmoothingScale<Dimension>::
zerothMoment() const {
  return mZerothMoment;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SPHSmoothingScale<Dimension>::
firstMoment() const {
  return mFirstMoment;
}

}
