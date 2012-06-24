namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// Access the polytropic constant.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
PolytropicEquationOfState<Dimension, Constants>::
polytropicConstant() const {
  return mPolytropicConstant;
}

//------------------------------------------------------------------------------
// Access the polytropic index.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
PolytropicEquationOfState<Dimension, Constants>::
polytropicIndex() const {
  return mPolytropicIndex;
}

//------------------------------------------------------------------------------
// Access gamma ( (n+1)/n ).
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
PolytropicEquationOfState<Dimension, Constants>::
gamma() const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Access the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
PolytropicEquationOfState<Dimension, Constants>::
molecularWeight() const {
  return mMolecularWeight;
}

//------------------------------------------------------------------------------
// Access the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
PolytropicEquationOfState<Dimension, Constants>::
externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension, typename Constants>
inline
void
PolytropicEquationOfState<Dimension, Constants>::
setExternalPressure(double P) {
  mExternalPressure = P;
}

}
}
