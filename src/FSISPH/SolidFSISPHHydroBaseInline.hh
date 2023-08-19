namespace Spheral {

//------------------------------------------------------------------------------
// set/get 
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
surfaceForceCoefficient(double x) {
  mSurfaceForceCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
surfaceForceCoefficient() const {
  return mSurfaceForceCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
densityStabilizationCoefficient(double x) {
  mDensityStabilizationCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
densityStabilizationCoefficient() const {
  return mDensityStabilizationCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
specificThermalEnergyDiffusionCoefficient(double x) {
  mSpecificThermalEnergyDiffusionCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
specificThermalEnergyDiffusionCoefficient() const {
  return mSpecificThermalEnergyDiffusionCoefficient;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
xsphCoefficient(double x) {
  mXSPHCoefficient = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
xsphCoefficient() const {
  return mXSPHCoefficient;
}


//------------------------------------------------------------------------------
// return our interface method
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
interfaceMethod(InterfaceMethod x) {
  mInterfaceMethod = x;
}
template<typename Dimension>
inline
InterfaceMethod
SolidFSISPHHydroBase<Dimension>::
interfaceMethod() const {
  return mInterfaceMethod;
}

//------------------------------------------------------------------------------
// return our kernel method
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
kernelAveragingMethod(KernelAveragingMethod x) {
  mKernelAveragingMethod = x;
}
template<typename Dimension>
inline
KernelAveragingMethod
SolidFSISPHHydroBase<Dimension>::
kernelAveragingMethod() const {
  return mKernelAveragingMethod;
}


//------------------------------------------------------------------------------
// switch to turn on density sum for different nodeLists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
applySelectSumDensity(bool x) {
  mApplySelectDensitySum = x;
}
template<typename Dimension>
inline
bool
SolidFSISPHHydroBase<Dimension>::
applySelectSumDensity() const {
  return mApplySelectDensitySum;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
sumDensityNodeLists(std::vector<int> x) {
  mSumDensityNodeLists = x;
}
template<typename Dimension>
inline
std::vector<int>
SolidFSISPHHydroBase<Dimension>::
sumDensityNodeLists() const {
  return mSumDensityNodeLists;
}

//------------------------------------------------------------------------------
// Return ref to our pair-wise energy derivs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
pairDepsDt() const {
  return mPairDepsDt;
}

//------------------------------------------------------------------------------
// Return the pressure gradient field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
DPDx() const {
  return mDPDx;
}

//------------------------------------------------------------------------------
// Return the specific thermal energy gradient field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
DepsDx() const {
  return mDepsDx;
}


//------------------------------------------------------------------------------
// our interface smoothness metric
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
rawPressure() const {
  return mRawPressure;
}

//------------------------------------------------------------------------------
// Return the interface normal field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
interfaceNormals() const {
  return mInterfaceNormals;
}

//------------------------------------------------------------------------------
// Return the Interface fraction
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceFraction() const {
  return mInterfaceFraction;
}

//------------------------------------------------------------------------------
// smoothness metric for mixing interfaces
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceSmoothness() const {
  return mInterfaceSmoothness;
}

//------------------------------------------------------------------------------
// next time step  Interface normal field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceNormals() const {
  return mNewInterfaceNormals;
}

//------------------------------------------------------------------------------
// next time step  Interface normal field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
smoothedInterfaceNormals() const {
  return mSmoothedInterfaceNormals;
}

//------------------------------------------------------------------------------
// next time step  Interface fraction
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceFraction() const {
  return mNewInterfaceFraction;
}

//------------------------------------------------------------------------------
// next time step smoothness metric for mixing interfaces
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceSmoothness() const {
  return mNewInterfaceSmoothness;
}


//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SolidFSISPHHydroBase<Dimension>::
DdeviatoricStressDt() const {
  return mDdeviatoricStressDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
bulkModulus() const {
  return mBulkModulus;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
shearModulus() const {
  return mShearModulus;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
yieldStrength() const {
  return mYieldStrength;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
plasticStrain0() const {
  return mPlasticStrain0;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SolidFSISPHHydroBase<Dimension>::
Hfield0() const {
  return mHfield0;
}


//------------------------------------------------------------------------------
// Ref to the slide surface obj
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SlideSurface<Dimension>&
SolidFSISPHHydroBase<Dimension>::slideSurface() const {
  return mSlideSurface;
}




//////////////////////////////////////////////////////////////////////////////////
// SPH HYDRO STUPH
//////////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------------
// Choose whether we want to sum for mass density, or integrate the continuity
// equation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
MassDensityType
SPHHydroBase<Dimension>::densityUpdate() const {
  return mDensityUpdate;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::
densityUpdate(MassDensityType type) {
  mDensityUpdate = type;
}

//------------------------------------------------------------------------------
// Choose how we want to update the H tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
HEvolutionType
SPHHydroBase<Dimension>::HEvolution() const {
  return mHEvolution;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::
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
SPHHydroBase<Dimension>::compatibleEnergyEvolution() const {
  return mCompatibleEnergyEvolution;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::compatibleEnergyEvolution(bool val) {
  mCompatibleEnergyEvolution = val;
}

//------------------------------------------------------------------------------
// Access the flag determining if we're evolving total or specific energy
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
SPHHydroBase<Dimension>::evolveTotalEnergy() const {
  return mEvolveTotalEnergy;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::evolveTotalEnergy(bool val) {
  mEvolveTotalEnergy = val;
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
SPHHydroBase<Dimension>::gradhCorrection(bool val) {
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
SPHHydroBase<Dimension>::XSPH(bool val) {
  mXSPH = val;
}

//------------------------------------------------------------------------------
// Access the flag controlling linear correct velocity gradient.
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
SPHHydroBase<Dimension>::correctVelocityGradient(bool val) {
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
SPHHydroBase<Dimension>::sumMassDensityOverAllNodeLists(bool val) {
  mSumMassDensityOverAllNodeLists = val;
}

//------------------------------------------------------------------------------
// Fraction of the centroidal filtering to apply.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
SPHHydroBase<Dimension>::filter() const {
  return mfilter;
}

template<typename Dimension>
inline
void
SPHHydroBase<Dimension>::filter(double val) {
  VERIFY(val >= 0.0 and val <= 1.0);
  mfilter = val;
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
epsilonTensile(typename Dimension::Scalar val) {
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
nTensile(typename Dimension::Scalar val) {
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
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
SPHHydroBase<Dimension>::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// Access the kernel used for artificial viscosity gradients.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
SPHHydroBase<Dimension>::
PiKernel() const {
  return mPiKernel;
}

//------------------------------------------------------------------------------
// The object defining how smoothing scales are evolved.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const SmoothingScaleBase<Dimension>&
SPHHydroBase<Dimension>::
smoothingScaleMethod() const {
  return mSmoothingScaleMethod;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
SPHHydroBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
pressure() const {
  return mPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
soundSpeed() const {
  return mSoundSpeed;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
omegaGradh() const {
  return mOmegaGradh;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
specificThermalEnergy0() const {
  return mSpecificThermalEnergy0;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
entropy() const {
  return mEntropy;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SPHHydroBase<Dimension>::
Hideal() const {
  return mHideal;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
maxViscousPressure() const {
  return mMaxViscousPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
effectiveViscousPressure() const {
  return mEffViscousPressure;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
massDensityCorrection() const {
  return mMassDensityCorrection;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
viscousWork() const {
  return mViscousWork;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
massDensitySum() const {
  return mMassDensitySum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
normalization() const {
  return mNormalization;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
weightedNeighborSum() const {
  return mWeightedNeighborSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SPHHydroBase<Dimension>::
massSecondMoment() const {
  return mMassSecondMoment;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
XSPHWeightSum() const {
  return mXSPHWeightSum;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SPHHydroBase<Dimension>::
XSPHDeltaV() const {
  return mXSPHDeltaV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SPHHydroBase<Dimension>::
M() const {
  return mM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SPHHydroBase<Dimension>::
localM() const {
  return mLocalM;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SPHHydroBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
SPHHydroBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
DmassDensityDt() const {
  return mDmassDensityDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
SPHHydroBase<Dimension>::
DspecificThermalEnergyDt() const {
  return mDspecificThermalEnergyDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::SymTensor>&
SPHHydroBase<Dimension>::
DHDt() const {
  return mDHDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SPHHydroBase<Dimension>::
DvDx() const {
  return mDvDx;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
SPHHydroBase<Dimension>::
internalDvDx() const {
  return mInternalDvDx;
}

template<typename Dimension>
inline
const std::vector<typename Dimension::Vector>&
SPHHydroBase<Dimension>::
pairAccelerations() const {
  return mPairAccelerations;
}

}
