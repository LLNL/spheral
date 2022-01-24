namespace Spheral {

//------------------------------------------------------------------------------
// Access the CFL safety criteria.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
GSPHHydroBase<Dimension>::cfl() const {
  return mCfl;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
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
GSPHHydroBase<Dimension>::useVelocityMagnitudeForDt() const {
  return mUseVelocityMagnitudeForDt;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
useVelocityMagnitudeForDt(bool x) {
  mUseVelocityMagnitudeForDt = x;
}

//------------------------------------------------------------------------------
// Return the master neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GSPHHydroBase<Dimension>::minMasterNeighbor() const {
  return mMinMasterNeighbor;
}

template<typename Dimension>
inline
int
GSPHHydroBase<Dimension>::maxMasterNeighbor() const {
  return mMaxMasterNeighbor;
}

template<typename Dimension>
inline
double
GSPHHydroBase<Dimension>::averageMasterNeighbor() const {
  return double(mSumMasterNeighbor)/(mNormMasterNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the coarse neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GSPHHydroBase<Dimension>::minCoarseNeighbor() const {
  return mMinCoarseNeighbor;
}

template<typename Dimension>
inline
int
GSPHHydroBase<Dimension>::maxCoarseNeighbor() const {
  return mMaxCoarseNeighbor;
}

template<typename Dimension>
inline
double
GSPHHydroBase<Dimension>::averageCoarseNeighbor() const {
  return double(mSumCoarseNeighbor)/(mNormCoarseNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the refine neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GSPHHydroBase<Dimension>::minRefineNeighbor() const {
  return mMinRefineNeighbor;
}

template<typename Dimension>
inline
int
GSPHHydroBase<Dimension>::maxRefineNeighbor() const {
  return mMaxRefineNeighbor;
}

template<typename Dimension>
inline
double
GSPHHydroBase<Dimension>::averageRefineNeighbor() const {
  return double(mSumRefineNeighbor)/(mNormRefineNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Return the actual neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GSPHHydroBase<Dimension>::minActualNeighbor() const {
  return mMinActualNeighbor;
}

template<typename Dimension>
inline
int
GSPHHydroBase<Dimension>::maxActualNeighbor() const {
  return mMaxActualNeighbor;
}

template<typename Dimension>
inline
double
GSPHHydroBase<Dimension>::averageActualNeighbor() const {
  return double(mSumActualNeighbor)/(mNormActualNeighbor + FLT_MIN);
}

//------------------------------------------------------------------------------
// Update the master neighboring statistics.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
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
GSPHHydroBase<Dimension>::
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
GSPHHydroBase<Dimension>::
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
GSPHHydroBase<Dimension>::
updateActualNeighborStats(int numNeighbor) const {
  if (numNeighbor > 0) {
    mMinActualNeighbor = std::min(mMinActualNeighbor, numNeighbor);
    mMaxActualNeighbor = std::max(mMaxActualNeighbor, numNeighbor);
    mSumActualNeighbor += numNeighbor;
    mNormActualNeighbor++;
  }
}

//------------------------------------------------------------------------------
// get out reimann solver obj
//------------------------------------------------------------------------------
template<typename Dimension>
inline
RiemannSolverBase<Dimension>&
GSPHHydroBase<Dimension>::riemannSolver() const {
  return mRiemannSolver;
}

//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MassDensityType
GSPHHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
densityUpdate(MassDensityType type) {
  mDensityUpdate = type;
}


//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
HEvolutionType
GSPHHydroBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
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
GSPHHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::compatibleEnergyEvolution(bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're evolving total or specific energy
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GSPHHydroBase<Dimension>::evolveTotalEnergy() const {
  return mEvolveTotalEnergy;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::evolveTotalEnergy(bool val) {
  mEvolveTotalEnergy = val;
}


//------------------------------------------------------------------------------
// Access the flag determining if we're using the XSPH algorithm.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GSPHHydroBase<Dimension>::XSPH() const {
  return mXSPH;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::XSPH(bool val) {
  mXSPH = val;
}

//------------------------------------------------------------------------------
// Access the flag controlling linear correct velocity gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GSPHHydroBase<Dimension>::correctVelocityGradient() const {
  return mCorrectVelocityGradient;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::correctVelocityGradient(bool val) {
  mCorrectVelocityGradient = val;
}


//------------------------------------------------------------------------------
// Parameter to determine the magnitude of the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
GSPHHydroBase<Dimension>::
epsilonTensile() const {
  return mEpsTensile;
}

template<typename Dimension>
void
GSPHHydroBase<Dimension>::
epsilonTensile(typename Dimension::Scalar val) {
  mEpsTensile = val;
}

//------------------------------------------------------------------------------
// Parameter to set the exponent used in the tensile small scale correction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
GSPHHydroBase<Dimension>::
nTensile() const {
  return mnTensile;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
nTensile(typename Dimension::Scalar val) {
  mnTensile = val;
}

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
GSPHHydroBase<Dimension>::
xmin() const {
  return mxmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
GSPHHydroBase<Dimension>::
xmax() const {
  return mxmax;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mxmin = x;
}

template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mxmax = x;
}

//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
GSPHHydroBase<Dimension>::
kernel() const {
  return mKernel;
}


//------------------------------------------------------------------------------
// The object defining how smoothing scales are evolved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const SmoothingScaleBase<Dimension>&
GSPHHydroBase<Dimension>::
smoothingScaleMethod() const {
  return mSmoothingScaleMethod;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
GSPHHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GSPHHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GSPHHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GSPHHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
GSPHHydroBase<Dimension>::
Hideal() const {
  return mHideal;
}


template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GSPHHydroBase<Dimension>::
normalization() const {
  return mNormalization;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GSPHHydroBase<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
GSPHHydroBase<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GSPHHydroBase<Dimension>::
XSPHWeightSum() const {
  return mXSPHWeightSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GSPHHydroBase<Dimension>::
XSPHDeltaV() const {
  return mXSPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GSPHHydroBase<Dimension>::
M() const {
  return mM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GSPHHydroBase<Dimension>::
localM() const {
  return mLocalM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GSPHHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GSPHHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GSPHHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
GSPHHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
GSPHHydroBase<Dimension>::
DHDt() const {
  return mDHDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GSPHHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GSPHHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GSPHHydroBase<Dimension>::
DpDx() const {
  return mDpDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
GSPHHydroBase<Dimension>::
DpDxRaw() const {
  return mDpDxRaw;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
GSPHHydroBase<Dimension>::
DvDxRaw() const {
  return mDvDxRaw;
}



// template<typename Dimension>
// inline
// const FieldList<Dimension, typename Dimension::Vector>&
// GSPHHydroBase<Dimension>::
// DrhoDx() const {
//   return mDrhoDx;
// }

template<typename Dimension>
inline
const std::vector<typename Dimension::Vector>&
GSPHHydroBase<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}


template<typename Dimension>
inline
void
GSPHHydroBase<Dimension>::
specificThermalEnergyDiffusionCoefficient(const typename Dimension::Scalar x) {
  mSpecificThermalEnergyDiffusionCoefficient = x;
}
template<typename Dimension>
inline
typename Dimension::Scalar
GSPHHydroBase<Dimension>::
specificThermalEnergyDiffusionCoefficient() const {
  return mSpecificThermalEnergyDiffusionCoefficient;
}

}
