namespace Spheral {
namespace SolidMaterial {

//------------------------------------------------------------------------------
// Access the coefficients. (n)
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
MurnahanEquationOfState<Dimension, Constants>::
n() const {
  return mn;
}

template<typename Dimension, typename Constants>
inline
void
MurnahanEquationOfState<Dimension, Constants>::
n(const double x) {
  mn = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (K)
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
MurnahanEquationOfState<Dimension, Constants>::
K() const {
  return mK;
}

template<typename Dimension, typename Constants>
inline
void
MurnahanEquationOfState<Dimension, Constants>::
K(const double x) {
  mK = x;
}

//------------------------------------------------------------------------------
// Access the atomic weight
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
MurnahanEquationOfState<Dimension, Constants>::
atomicWeight() const {
  return mAtomicWeight;
}

template<typename Dimension, typename Constants>
inline
void
MurnahanEquationOfState<Dimension, Constants>::
atomicWeight(const double x) {
  mAtomicWeight = x;
}

//------------------------------------------------------------------------------
// Get and set the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
MurnahanEquationOfState<Dimension, Constants>::
externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension, typename Constants>
inline
void
MurnahanEquationOfState<Dimension, Constants>::
externalPressure(const double val) {
  mExternalPressure = val;
}

}
}
