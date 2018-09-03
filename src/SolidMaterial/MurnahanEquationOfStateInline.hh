namespace Spheral {

//------------------------------------------------------------------------------
// Access the coefficients. (n)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
MurnahanEquationOfState<Dimension>::
n() const {
  return mn;
}

template<typename Dimension>
inline
void
MurnahanEquationOfState<Dimension>::
n(const double x) {
  mn = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (K)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
MurnahanEquationOfState<Dimension>::
K() const {
  return mK;
}

template<typename Dimension>
inline
void
MurnahanEquationOfState<Dimension>::
K(const double x) {
  mK = x;
}

//------------------------------------------------------------------------------
// Access the atomic weight
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
MurnahanEquationOfState<Dimension>::
atomicWeight() const {
  return mAtomicWeight;
}

template<typename Dimension>
inline
void
MurnahanEquationOfState<Dimension>::
atomicWeight(const double x) {
  mAtomicWeight = x;
}

//------------------------------------------------------------------------------
// Get and set the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
MurnahanEquationOfState<Dimension>::
externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension>
inline
void
MurnahanEquationOfState<Dimension>::
externalPressure(const double val) {
  mExternalPressure = val;
}

}
