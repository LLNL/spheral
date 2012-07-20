namespace Spheral {
namespace SolidMaterial {

//------------------------------------------------------------------------------
// Get and set C0.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
GruneisenEquationOfState<Dimension, Constants>::C0() const {
  return mC0;
}

template<typename Dimension, typename Constants>
inline
void
GruneisenEquationOfState<Dimension, Constants>::C0(const double C0) {
  mC0 = C0;
}

//------------------------------------------------------------------------------
// Get and set S1.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
GruneisenEquationOfState<Dimension, Constants>::S1() const {
  return mS1;
}

template<typename Dimension, typename Constants>
inline
void
GruneisenEquationOfState<Dimension, Constants>::S1(const double S1) {
  mS1 = S1;
}

//------------------------------------------------------------------------------
// Get and set S2.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
GruneisenEquationOfState<Dimension, Constants>::S2() const {
  return mS2;
}

template<typename Dimension, typename Constants>
inline
void
GruneisenEquationOfState<Dimension, Constants>::S2(const double S2) {
  mS2 = S2;
}

//------------------------------------------------------------------------------
// Get and set S3.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
GruneisenEquationOfState<Dimension, Constants>::S3() const {
  return mS3;
}

template<typename Dimension, typename Constants>
inline
void
GruneisenEquationOfState<Dimension, Constants>::S3(const double S3) {
  mS3 = S3;
}

//------------------------------------------------------------------------------
// Get and set gamma0.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
GruneisenEquationOfState<Dimension, Constants>::gamma0() const {
  return mgamma0;
}

template<typename Dimension, typename Constants>
inline
void
GruneisenEquationOfState<Dimension, Constants>::gamma0(const double gamma0) {
  mgamma0 = gamma0;
}

//------------------------------------------------------------------------------
// Get and set b.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
GruneisenEquationOfState<Dimension, Constants>::b() const {
  return mb;
}

template<typename Dimension, typename Constants>
inline
void
GruneisenEquationOfState<Dimension, Constants>::b(const double b) {
  mb = b;
}

//------------------------------------------------------------------------------
// Get and set the atomic weight.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
GruneisenEquationOfState<Dimension, Constants>::atomicWeight() const {
  return mAtomicWeight;
}

template<typename Dimension, typename Constants>
inline
void
GruneisenEquationOfState<Dimension, Constants>::atomicWeight(const double val) {
  mAtomicWeight = val;
}

//------------------------------------------------------------------------------
// Get Cv.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
GruneisenEquationOfState<Dimension, Constants>::Cv() const {
  return mCv;
}

//------------------------------------------------------------------------------
// Get and set the external pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
inline
double
GruneisenEquationOfState<Dimension, Constants>::externalPressure() const {
  return mExternalPressure;
}

template<typename Dimension, typename Constants>
inline
void
GruneisenEquationOfState<Dimension, Constants>::externalPressure(const double val) {
  mExternalPressure = val;
}

}
}
