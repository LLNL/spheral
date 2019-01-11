namespace Spheral {

//------------------------------------------------------------------------------
// Access the polytropic constant.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
IsothermalEquationOfState<Dimension>::
K() const {
  return mK;
}

//------------------------------------------------------------------------------
// Access the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
IsothermalEquationOfState<Dimension>::
molecularWeight() const {
  return mMolecularWeight;
}

//------------------------------------------------------------------------------
// Access the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
IsothermalEquationOfState<Dimension>::
externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension>
inline
void
IsothermalEquationOfState<Dimension>::
setExternalPressure(double P) {
  mExternalPressure = P;
}

}
