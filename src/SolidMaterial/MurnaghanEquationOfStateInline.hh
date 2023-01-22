namespace Spheral {

//------------------------------------------------------------------------------
// Access the coefficients. (n)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
MurnaghanEquationOfState<Dimension>::
n() const {
  return mn;
}

template<typename Dimension>
inline
void
MurnaghanEquationOfState<Dimension>::
n(double x) {
  mn = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (K)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
MurnaghanEquationOfState<Dimension>::
K() const {
  return mK;
}

template<typename Dimension>
inline
void
MurnaghanEquationOfState<Dimension>::
K(double x) {
  mK = x;
}

//------------------------------------------------------------------------------
// Access the atomic weight
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
MurnaghanEquationOfState<Dimension>::
atomicWeight() const {
  return mAtomicWeight;
}

template<typename Dimension>
inline
void
MurnaghanEquationOfState<Dimension>::
atomicWeight(double x) {
  mAtomicWeight = x;
}

}
