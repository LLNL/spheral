namespace Spheral {

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
ASPHSmoothingScale<Dimension>::
WT() const {
  return mWT;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
ASPHSmoothingScale<Dimension>::
zerothMoment() const {
  return mZerothMoment;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
ASPHSmoothingScale<Dimension>::
firstMoment() const {
  return mFirstMoment;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
ASPHSmoothingScale<Dimension>::
secondMoment() const {
  return mSecondMoment;
}

}
