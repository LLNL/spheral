namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// Access the polytropic constant.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
IsothermalEquationOfState<Dimension, Constants>::
K() const {
  return mK;
}

//------------------------------------------------------------------------------
// Access the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
IsothermalEquationOfState<Dimension, Constants>::
molecularWeight() const {
  return mMolecularWeight;
}

//------------------------------------------------------------------------------
// Access the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
IsothermalEquationOfState<Dimension, Constants>::
externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension, typename Constants>
inline
void
IsothermalEquationOfState<Dimension, Constants>::
setExternalPressure(double P) {
  mExternalPressure = P;
}

}
}
