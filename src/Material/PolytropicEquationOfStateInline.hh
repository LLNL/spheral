namespace Spheral {

//------------------------------------------------------------------------------
// Access the polytropic constant.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PolytropicEquationOfState<Dimension>::
polytropicConstant() const {
  return mPolytropicConstant;
}

//------------------------------------------------------------------------------
// Access the polytropic index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PolytropicEquationOfState<Dimension>::
polytropicIndex() const {
  return mPolytropicIndex;
}

//------------------------------------------------------------------------------
// Access gamma ( (n+1)/n ).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PolytropicEquationOfState<Dimension>::
gamma() const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Access the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PolytropicEquationOfState<Dimension>::
molecularWeight() const {
  return mMolecularWeight;
}

//------------------------------------------------------------------------------
// Access the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
PolytropicEquationOfState<Dimension>::
externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension>
inline
void
PolytropicEquationOfState<Dimension>::
setExternalPressure(double P) {
  mExternalPressure = P;
}

}
