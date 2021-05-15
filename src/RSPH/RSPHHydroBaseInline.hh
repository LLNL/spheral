namespace Spheral {

////////////////////////////////////////////////////////////////////////////////
// Generic Hydro Methods
////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
// Access the CFL safety criteria.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
RSPHHydroBase<Dimension>::cfl() const {
  return mCfl;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
cfl(typename Dimension::Scalar cfl) {
  mCfl = cfl;
}

//------------------------------------------------------------------------------
// Flag to set whether or not to use the magnitude of the velocity to set the
// timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RSPHHydroBase<Dimension>::useVelocityMagnitudeForDt() const {
  return mUseVelocityMagnitudeForDt;
}

template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
useVelocityMagnitudeForDt(bool x) {
  mUseVelocityMagnitudeForDt = x;
}

//------------------------------------------------------------------------------
// Return the master neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RSPHHydroBase<Dimension>::minMasterNeighbor() const {
  return mMinMasterNeighbor;
}

template<typename Dimension>
inline
int
RSPHHydroBase<Dimension>::maxMasterNeighbor() const {
  return mMaxMasterNeighbor;
}

template<typename Dimension>
inline
double
RSPHHydroBase<Dimension>::averageMasterNeighbor() const {
  return double(mSumMasterNeighbor)/(mNormMasterNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the coarse neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RSPHHydroBase<Dimension>::minCoarseNeighbor() const {
  return mMinCoarseNeighbor;
}

template<typename Dimension>
inline
int
RSPHHydroBase<Dimension>::maxCoarseNeighbor() const {
  return mMaxCoarseNeighbor;
}

template<typename Dimension>
inline
double
RSPHHydroBase<Dimension>::averageCoarseNeighbor() const {
  return double(mSumCoarseNeighbor)/(mNormCoarseNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the refine neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RSPHHydroBase<Dimension>::minRefineNeighbor() const {
  return mMinRefineNeighbor;
}

template<typename Dimension>
inline
int
RSPHHydroBase<Dimension>::maxRefineNeighbor() const {
  return mMaxRefineNeighbor;
}

template<typename Dimension>
inline
double
RSPHHydroBase<Dimension>::averageRefineNeighbor() const {
  return double(mSumRefineNeighbor)/(mNormRefineNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the actual neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RSPHHydroBase<Dimension>::minActualNeighbor() const {
  return mMinActualNeighbor;
}

template<typename Dimension>
inline
int
RSPHHydroBase<Dimension>::maxActualNeighbor() const {
  return mMaxActualNeighbor;
}

template<typename Dimension>
inline
double
RSPHHydroBase<Dimension>::averageActualNeighbor() const {
  return double(mSumActualNeighbor)/(mNormActualNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Update the master neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
updateMasterNeighborStats(int numMaster) const {
  if (numMaster > 0) {
    mMinMasterNeighbor = std::min(mMinMasterNeighbor, numMaster);
    mMaxMasterNeighbor = std::max(mMaxMasterNeighbor, numMaster);
    mSumMasterNeighbor += numMaster;
    mNormMasterNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the coarse neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
updateCoarseNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinCoarseNeighbor = std::min(mMinCoarseNeighbor, numNeighbor);
    mMaxCoarseNeighbor = std::max(mMaxCoarseNeighbor, numNeighbor);
    mSumCoarseNeighbor += numNeighbor;
    mNormCoarseNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the refine neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
updateRefineNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinRefineNeighbor = std::min(mMinRefineNeighbor, numNeighbor);
    mMaxRefineNeighbor = std::max(mMaxRefineNeighbor, numNeighbor);
    mSumRefineNeighbor += numNeighbor;
    mNormRefineNeighbor++;
  }
}

//------------------------------------------------------------------------------
// Update the actual neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
RSPHHydroBase<Dimension>::
updateActualNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinActualNeighbor = std::min(mMinActualNeighbor, numNeighbor);
    mMaxActualNeighbor = std::max(mMaxActualNeighbor, numNeighbor);
    mSumActualNeighbor += numNeighbor;
    mNormActualNeighbor++;
  }
}










////////////////////////////////////////////////////////////////////////////////
// SPH Hydro Methods
////////////////////////////////////////////////////////////////////////////////
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
specificThermalEnergy0() const {
  return mSpecificThermalEnergy0;
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
