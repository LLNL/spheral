namespace Spheral {

//------------------------------------------------------------------------------
// a1
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::a1() const {
  return mA1;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::a1(const double val) {
  mA1 = val;
}

//------------------------------------------------------------------------------
// a2pos
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::a2pos() const {
  return mA2pos;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::a2pos(const double val) {
  mA2pos = val;
}

//------------------------------------------------------------------------------
// a2neg
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::a2neg() const {
  return mA2neg;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::a2neg(const double val) {
  mA2neg = val;
}

//------------------------------------------------------------------------------
// b0
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::b0() const {
  return mB0;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::b0(const double val) {
  mB0 = val;
}

//------------------------------------------------------------------------------
// b1
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::b1() const {
  return mB1;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::b1(const double val) {
  mB1 = val;
}

//------------------------------------------------------------------------------
// b2pos
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::b2pos() const {
  return mB2pos;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::b2pos(const double val) {
  mB2pos = val;
}

//------------------------------------------------------------------------------
// b2neg
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::b2neg() const {
  return mB2neg;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::b2neg(const double val) {
  mB2neg = val;
}

//------------------------------------------------------------------------------
// c0
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::c0() const {
  return mC0;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::c0(const double val) {
  mC0 = val;
}

//------------------------------------------------------------------------------
// c1
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::c1() const {
  return mC1;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::c1(const double val) {
  mC1 = val;
}

//------------------------------------------------------------------------------
// c2pos
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::c2pos() const {
  return mC2pos;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::c2pos(const double val) {
  mC2pos = val;
}

//------------------------------------------------------------------------------
// c2neg
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::c2neg() const {
  return mC2neg;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::c2neg(const double val) {
  mC2neg = val;
}

//------------------------------------------------------------------------------
// E0
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::E0() const {
  return mE0;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::E0(const double val) {
  mE0 = val;
}

//------------------------------------------------------------------------------
// atomicWeight
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::atomicWeight() const {
  return mAtomicWeight;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::atomicWeight(const double val) {
  mAtomicWeight = val;
}

//------------------------------------------------------------------------------
// Get Cv.
//------------------------------------------------------------------------------
template<typename Dimension>
double
OsborneEquationOfState<Dimension>::Cv() const {
  return mCv;
}

//------------------------------------------------------------------------------
// Get and set the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
OsborneEquationOfState<Dimension>::externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension>
inline
void
OsborneEquationOfState<Dimension>::externalPressure(const double val) {
  mExternalPressure = val;
}

}
