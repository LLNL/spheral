namespace Spheral {

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (A0)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
LinearPolynomialEquationOfState<Dimension>::
a0() const {
  return mA0;
}

template<typename Dimension>
inline
void
LinearPolynomialEquationOfState<Dimension>::
a0(double x) {
  mA0 = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (A1)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
LinearPolynomialEquationOfState<Dimension>::
a1() const {
  return mA1;
}

template<typename Dimension>
inline
void
LinearPolynomialEquationOfState<Dimension>::
a1(double x) {
  mA1 = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (A2)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
LinearPolynomialEquationOfState<Dimension>::
a2() const {
  return mA2;
}

template<typename Dimension>
inline
void
LinearPolynomialEquationOfState<Dimension>::
a2(double x) {
  mA2 = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (A3)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
LinearPolynomialEquationOfState<Dimension>::
a3() const {
  return mA3;
}

template<typename Dimension>
inline
void
LinearPolynomialEquationOfState<Dimension>::
a3(double x) {
  mA3 = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (B0)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
LinearPolynomialEquationOfState<Dimension>::
b0() const {
  return mB0;
}

template<typename Dimension>
inline
void
LinearPolynomialEquationOfState<Dimension>::
b0(double x) {
  mB0 = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (B1)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
LinearPolynomialEquationOfState<Dimension>::
b1() const {
  return mB1;
}

template<typename Dimension>
inline
void
LinearPolynomialEquationOfState<Dimension>::
b1(double x) {
  mB1 = x;
}

//------------------------------------------------------------------------------
// Access the polynomial coefficients. (B2)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
LinearPolynomialEquationOfState<Dimension>::
b2() const {
  return mB2;
}

template<typename Dimension>
inline
void
LinearPolynomialEquationOfState<Dimension>::
b2(double x) {
  mB2 = x;
}

//------------------------------------------------------------------------------
// Access the atomic weight
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
LinearPolynomialEquationOfState<Dimension>::
atomicWeight() const {
  return mAtomicWeight;
}

template<typename Dimension>
inline
void
LinearPolynomialEquationOfState<Dimension>::
atomicWeight(double x) {
  mAtomicWeight = x;
}

}
