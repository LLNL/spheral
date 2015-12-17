namespace Spheral {
namespace SPHSpace {

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
PSPHHydroBase<Dimension>::
gamma() const {
  return mGamma;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
PSPHHydroBase<Dimension>::
PSPHpbar() const {
  return mPSPHpbar;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
PSPHHydroBase<Dimension>::
PSPHcorrection() const {
  return mPSPHcorrection;
}

}
}
