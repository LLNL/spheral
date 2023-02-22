namespace Spheral {

//------------------------------------------------------------------------------
// set/get 
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
planeStrain(bool x) {
  mPlaneStrain = x;
}
template<typename Dimension>
inline
bool
SolidFSISPHHydroBase<Dimension>::
planeStrain() const {
  return mPlaneStrain;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
interfacePmin(double x) {
  mInterfacePmin = x;
}
template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
interfacePmin() const {
  return mInterfacePmin;
}

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
// next time step interface flags
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  int>&
SolidFSISPHHydroBase<Dimension>::
interfaceFlags() const {
  return mInterfaceFlags;
}

//------------------------------------------------------------------------------
// Return the interface area vectors field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
interfaceAreaVectors() const {
  return mInterfaceAreaVectors;
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
// next time step interface flags
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  int>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceFlags() const {
  return mNewInterfaceFlags;
}

//------------------------------------------------------------------------------
// next time step interface area vectors field list ref
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Vector>&
SolidFSISPHHydroBase<Dimension>::
newInterfaceAreaVectors() const {
  return mNewInterfaceAreaVectors;
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
// normalization for our smoothness calculation
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceSmoothnessNormalization() const {
  return mInterfaceSmoothnessNormalization;
}

//------------------------------------------------------------------------------
// normalization for only same material nodes
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceFraction() const {
  return mInterfaceFraction;
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
// next time step smoothness metric for mixing interfaces
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension,  typename Dimension::Scalar>&
SolidFSISPHHydroBase<Dimension>::
interfaceAngles() const {
  return mInterfaceAngles;
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

}
