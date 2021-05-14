namespace Spheral {



//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
HEvolutionType
RSPHHydroBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
HEvolution(HEvolutionType type) {
  mHEvolution = type;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the compatible energy evolution 
// algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RSPHHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::compatibleEnergyEvolution(bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're evolving total or specific energy
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RSPHHydroBase<Dimension>::evolveTotalEnergy() const {
  return mEvolveTotalEnergy;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::evolveTotalEnergy(bool val) {
  mEvolveTotalEnergy = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the grad h correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RSPHHydroBase<Dimension>::gradhCorrection() const {
  return mGradhCorrection;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::gradhCorrection(bool val) {
  mGradhCorrection = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RSPHHydroBase<Dimension>::XSPH() const {
  return mXSPH;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::XSPH(bool val) {
  mXSPH = val;
}

//------------------------------------------------------------------------------
// Access the flag controlling linear correct velocity gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RSPHHydroBase<Dimension>::correctVelocityGradient() const {
  return mCorrectVelocityGradient;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::correctVelocityGradient(bool val) {
  mCorrectVelocityGradient = val;
}


//------------------------------------------------------------------------------
// Fraction of the centroidal filtering to apply.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
RSPHHydroBase<Dimension>::filter() const {
  return mfilter;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::filter(double val) {
  VERIFY(val >= 0.0 and val <= 1.0);
  mfilter = val;
}

//------------------------------------------------------------------------------
// Parameter to determine the magnitude of the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
RSPHHydroBase<Dimension>::
epsilonTensile() const {
  return mEpsTensile;
}

template<typename Dimension>
void
RSPHHydroBase<Dimension>::
epsilonTensile(typename Dimension::Scalar val) {
  mEpsTensile = val;
}

//------------------------------------------------------------------------------
// Parameter to set the exponent used in the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
RSPHHydroBase<Dimension>::
nTensile() const {
  return mnTensile;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
nTensile(typename Dimension::Scalar val) {
  mnTensile = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
RSPHHydroBase<Dimension>::
xmin() const {
  return mxmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
RSPHHydroBase<Dimension>::
xmax() const {
  return mxmax;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mxmin = x;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mxmax = x;
}

//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
RSPHHydroBase<Dimension>::
kernel() const {
  return mKernel;
}


//------------------------------------------------------------------------------
// The object defining how smoothing scales are evolved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const SmoothingScaleBase<Dimension>&
RSPHHydroBase<Dimension>::
smoothingScaleMethod() const {
  return mSmoothingScaleMethod;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
RSPHHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
omegaGradh() const {
  return mOmegaGradh;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
specificThermalEnergy0() const {
  return mSpecificThermalEnergy0;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
entropy() const {
  return mEntropy;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
RSPHHydroBase<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
maxViscousPressure() const {
  return mMaxViscousPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
effectiveViscousPressure() const {
  return mEffViscousPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
massDensityCorrection() const {
  return mMassDensityCorrection;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
viscousWork() const {
  return mViscousWork;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
massDensitySum() const {
  return mMassDensitySum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
normalization() const {
  return mNormalization;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
RSPHHydroBase<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
XSPHWeightSum() const {
  return mXSPHWeightSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
RSPHHydroBase<Dimension>::
XSPHDeltaV() const {
  return mXSPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
RSPHHydroBase<Dimension>::
M() const {
  return mM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
RSPHHydroBase<Dimension>::
localM() const {
  return mLocalM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
RSPHHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
RSPHHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
RSPHHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
RSPHHydroBase<Dimension>::
DHDt() const {
  return mDHDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
RSPHHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
RSPHHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const std::vector<typename Dimension::Vector>&
RSPHHydroBase<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}

}
