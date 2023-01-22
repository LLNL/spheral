namespace Spheral {

//------------------------------------------------------------------------------
// Get and set C0.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GruneisenEquationOfState<Dimension>::C0() const {
  return mC0;
}

template<typename Dimension>
inline
void
GruneisenEquationOfState<Dimension>::C0(double C0) {
  mC0 = C0;
}

//------------------------------------------------------------------------------
// Get and set S1.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GruneisenEquationOfState<Dimension>::S1() const {
  return mS1;
}

template<typename Dimension>
inline
void
GruneisenEquationOfState<Dimension>::S1(double S1) {
  mS1 = S1;
}

//------------------------------------------------------------------------------
// Get and set S2.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GruneisenEquationOfState<Dimension>::S2() const {
  return mS2;
}

template<typename Dimension>
inline
void
GruneisenEquationOfState<Dimension>::S2(double S2) {
  mS2 = S2;
}

//------------------------------------------------------------------------------
// Get and set S3.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GruneisenEquationOfState<Dimension>::S3() const {
  return mS3;
}

template<typename Dimension>
inline
void
GruneisenEquationOfState<Dimension>::S3(double S3) {
  mS3 = S3;
}

//------------------------------------------------------------------------------
// Get and set gamma0.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GruneisenEquationOfState<Dimension>::gamma0() const {
  return mgamma0;
}

template<typename Dimension>
inline
void
GruneisenEquationOfState<Dimension>::gamma0(double gamma0) {
  mgamma0 = gamma0;
}

//------------------------------------------------------------------------------
// Get and set b.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GruneisenEquationOfState<Dimension>::b() const {
  return mb;
}

template<typename Dimension>
inline
void
GruneisenEquationOfState<Dimension>::b(double b) {
  mb = b;
}

//------------------------------------------------------------------------------
// Get and set the atomic weight.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GruneisenEquationOfState<Dimension>::atomicWeight() const {
  return mAtomicWeight;
}

template<typename Dimension>
inline
void
GruneisenEquationOfState<Dimension>::atomicWeight(double val) {
  mAtomicWeight = val;
}

//------------------------------------------------------------------------------
// Get and set the energy multiplier.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GruneisenEquationOfState<Dimension>::energyMultiplier() const {
  return mEnergyMultiplier;
}

template<typename Dimension>
inline
void
GruneisenEquationOfState<Dimension>::energyMultiplier(double val) {
  mEnergyMultiplier = val;
}

//------------------------------------------------------------------------------
// Get Cv.
//------------------------------------------------------------------------------
template<typename Dimension>
double
GruneisenEquationOfState<Dimension>::Cv() const {
  return mCv;
}

}
