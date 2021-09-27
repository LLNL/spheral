namespace Spheral {

//------------------------------------------------------------------------------
// Get and set gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
StiffenedGas<Dimension>::gamma() const {
  return mGamma;
}

template<typename Dimension>
inline
void
StiffenedGas<Dimension>::gamma(typename Dimension::Scalar gam) {
  mGamma = gam;
  mGamma1 = mGamma - 1.0;
}

//------------------------------------------------------------------------------
// Get and set the specific Heat.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
StiffenedGas<Dimension>::specificHeat() const {
  return mCv;
}

template<typename Dimension>
inline
void
StiffenedGas<Dimension>::specificHeat(typename Dimension::Scalar Cv) {
  mCv = Cv;
}


//------------------------------------------------------------------------------
// Get / Set referencePressure
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
StiffenedGas<Dimension>::referencePressure() const {
  return mP0;
}

template<typename Dimension>
inline
void
StiffenedGas<Dimension>::referencePressure(typename Dimension::Scalar P0) {
  mP0 = P0;
}

}
