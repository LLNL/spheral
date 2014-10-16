namespace Spheral {
namespace SPHSpace {

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
PhysicsSpace::MassDensityType
SPHHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::
densityUpdate(const PhysicsSpace::MassDensityType type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
PhysicsSpace::HEvolutionType
SPHHydroBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::
HEvolution(const PhysicsSpace::HEvolutionType type) {
  mHEvolution = type;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the compatible energy evolution 
// algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SPHHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::compatibleEnergyEvolution(const bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the grad h correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SPHHydroBase<Dimension>::gradhCorrection() const {
  return mGradhCorrection;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::gradhCorrection(const bool val) {
  mGradhCorrection = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SPHHydroBase<Dimension>::XSPH() const {
  return mXSPH;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::XSPH(const bool val) {
  mXSPH = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SPHHydroBase<Dimension>::correctVelocityGradient() const {
  return mCorrectVelocityGradient;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::correctVelocityGradient(const bool val) {
  mCorrectVelocityGradient = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if the sum mass density definition goes over
// neighbor NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SPHHydroBase<Dimension>::sumMassDensityOverAllNodeLists() const {
  return mSumMassDensityOverAllNodeLists;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::sumMassDensityOverAllNodeLists(const bool val) {
  mSumMassDensityOverAllNodeLists = val;
}

//------------------------------------------------------------------------------
// Parameter to determine the magnitude of the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SPHHydroBase<Dimension>::
epsilonTensile() const {
  return mEpsTensile;
}

template<typename Dimension>
void
SPHHydroBase<Dimension>::
epsilonTensile(const typename Dimension::Scalar val) {
  mEpsTensile = val;
}

//------------------------------------------------------------------------------
// Parameter to set the exponent used in the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SPHHydroBase<Dimension>::
nTensile() const {
  return mnTensile;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::
nTensile(const typename Dimension::Scalar val) {
  mnTensile = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
SPHHydroBase<Dimension>::
xmin() const {
  return mxmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SPHHydroBase<Dimension>::
xmax() const {
  return mxmax;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mxmin = x;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mxmax = x;
}

//------------------------------------------------------------------------------
// The object defining how smoothing scales are evolved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeSpace::SmoothingScaleBase<Dimension>&
SPHHydroBase<Dimension>::
smoothingScaleMethod() const {
  return mSmoothingScaleMethod;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, int>&
SPHHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
omegaGradh() const {
  return mOmegaGradh;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
specificThermalEnergy0() const {
  return mSpecificThermalEnergy0;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SPHHydroBase<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
maxViscousPressure() const {
  return mMaxViscousPressure;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
massDensitySum() const {
  return mMassDensitySum;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
normalization() const {
  return mNormalization;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SPHHydroBase<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
XSPHWeightSum() const {
  return mXSPHWeightSum;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SPHHydroBase<Dimension>::
XSPHDeltaV() const {
  return mXSPHDeltaV;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SPHHydroBase<Dimension>::
M() const {
  return mM;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SPHHydroBase<Dimension>::
localM() const {
  return mLocalM;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SPHHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Vector>&
SPHHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>&
SPHHydroBase<Dimension>::
DHDt() const {
  return mDHDt;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SPHHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
SPHHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, std::vector<typename Dimension::Vector> >&
SPHHydroBase<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}

}
}
