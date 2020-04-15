namespace Spheral {

//------------------------------------------------------------------------------
// Access the polytropic constant.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
K() const {
  return mK;
}

//------------------------------------------------------------------------------
// Access the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
molecularWeight() const {
  return mMolecularWeight;
}

//------------------------------------------------------------------------------
// Access the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension>
inline
void
IsothermalEquationOfState<Dimension>::
setExternalPressure(typename Dimension::Scalar P) {
  mExternalPressure = P;
}

}
