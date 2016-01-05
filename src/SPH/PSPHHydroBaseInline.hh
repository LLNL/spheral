namespace Spheral {
namespace SPHSpace {

//------------------------------------------------------------------------------
// Access the flag determining if we're applying Hopkins 2014 artificial
// conductivity.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
PSPHHydroBase<Dimension>::HopkinsConductivity() const {
  return mHopkinsConductivity;
}

template<typename Dimension>
inline
void
PSPHHydroBase<Dimension>::HopkinsConductivity(const bool val) {
  mHopkinsConductivity = val;
}

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
PSPHcorrection() const {
  return mPSPHcorrection;
}

}
}
