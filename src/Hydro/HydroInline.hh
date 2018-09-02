namespace Spheral {

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MassDensityType
Hydro<Dimension>::sumForMassDensity() const {
  return mSumForMassDensity;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::
sumForMassDensity(const MassDensityType type) {
  mSumForMassDensity = type;
}

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
HEvolutionType
Hydro<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::
HEvolution(const HEvolutionType type) {
  mHEvolution = type;
}

//------------------------------------------------------------------------------
// Access the minimum allowed smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
Hydro<Dimension>::hmin() const {
  CHECK(mhmin > 0.0);
  return mhmin;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::
hmin(const typename Dimension::Scalar val) {
  CHECK(val > 0.0);
  mhmin = val;
}

//------------------------------------------------------------------------------
// Access the maximum allowed smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
Hydro<Dimension>::hmax() const {
  CHECK(mhmax > 0.0);
  return mhmax;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::
hmax(const typename Dimension::Scalar val) {
  CHECK(val > 0.0);
  mhmax = val;
}

//------------------------------------------------------------------------------
// Access the minimum allowed ratio of the smoothing scales in the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
Hydro<Dimension>::hratiomin() const {
  CHECK(mhratiomin >= 0.0);
  return mhratiomin;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::
hratiomin(typename Dimension::Scalar val) {
  CHECK(val >= 0.0);
  mhratiomin = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the compatible energy evolution 
// algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Hydro<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::compatibleEnergyEvolution(const bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the grad h correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Hydro<Dimension>::gradhCorrection() const {
  return mGradhCorrection;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::gradhCorrection(const bool val) {
  mGradhCorrection = val;
}

//------------------------------------------------------------------------------
// Post iterate h number of cycles between firing.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
Hydro<Dimension>::postIterateHCycle() const {
  return mPostIterateHCycle;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::postIterateHCycle(const int val) {
  mPostIterateHCycle = val;
}

//------------------------------------------------------------------------------
// Post iterate h max iterations.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
Hydro<Dimension>::postIterateHMaxIterations() const {
  return mPostIterateHMaxIterations;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::postIterateHMaxIterations(const int val) {
  mPostIterateHMaxIterations = val;
}

//------------------------------------------------------------------------------
// Post iterate h tolerance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
Hydro<Dimension>::postIterateHtolerance() const {
  return mPostIterateHtolerance;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::postIterateHtolerance(const double val) {
  mPostIterateHtolerance = val;
}

//------------------------------------------------------------------------------
// Post iterate h n perh h.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
Hydro<Dimension>::postIterateHnPerh() const {
  return mPostIterateHnPerh;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::postIterateHnPerh(const double val) {
  mPostIterateHnPerh = val;
}

//------------------------------------------------------------------------------
// Post iterate h spherical start.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Hydro<Dimension>::postIterateHsphericalStart() const {
  return mPostIterateHsphericalStart;
}

template<typename Dimension>
inline
void
Hydro<Dimension>::postIterateHsphericalStart(const bool val) {
  mPostIterateHsphericalStart = val;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldList<Dimension, typename Dimension::SymTensor>
Hydro<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
FieldList<Dimension, int>
Hydro<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
FieldList<Dimension, typename Dimension::Scalar>
Hydro<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
FieldList<Dimension, typename Dimension::Scalar>
Hydro<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
FieldList<Dimension, typename Dimension::Scalar>
Hydro<Dimension>::
positionWeight() const {
  return mPositionWeight;
}

template<typename Dimension>
inline
FieldList<Dimension, std::vector<typename Dimension::Vector> >
Hydro<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}

template<typename Dimension>
inline
FieldList<Dimension, std::vector<typename Dimension::Vector> >
Hydro<Dimension>::
QpairAccelerations() const {
  return mQpairAccelerations;
}

template<typename Dimension>
inline
FieldList<Dimension, typename Dimension::Scalar>
Hydro<Dimension>::
specificThermalEnergy0() const {
  return mSpecificThermalEnergy0;
}

}
