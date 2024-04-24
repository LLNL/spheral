namespace Spheral {

//------------------------------------------------------------------------------
// Given the volume and target nperh, compute an effective target hmax
//------------------------------------------------------------------------------
// 1D
template<>
inline
typename Dim<1>::Scalar
SmoothingScaleBase<Dim<1>>::hmax(const Dim<1>::Scalar Vi,
                                 const Dim<1>::Scalar nperh) const {
  return 0.5*nperh*Vi;
}

// 2D
template<>
inline
typename Dim<1>::Scalar
SmoothingScaleBase<Dim<2>>::hmax(const Dim<2>::Scalar Vi,
                                 const Dim<2>::Scalar nperh) const {
  return nperh*std::sqrt(Vi/M_PI);
}

// 3D
template<>
inline
typename Dim<1>::Scalar
SmoothingScaleBase<Dim<3>>::hmax(const Dim<3>::Scalar Vi,
                                 const Dim<3>::Scalar nperh) const {
  return nperh*pow(0.75*Vi/M_PI, 1.0/3.0);
}

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
HEvolutionType
SmoothingScaleBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
SmoothingScaleBase<Dimension>::
HEvolution(HEvolutionType type) {
  mHEvolution = type;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SmoothingScaleBase<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SmoothingScaleBase<Dimension>::
DHDt() const {
  return mDHDt;
}

}
