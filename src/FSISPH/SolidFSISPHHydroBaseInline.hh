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
// return our slide surface method
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
slideSurfaceMethod(SlideSurfaceMethod x) {
  mSlideSurfaceMethod = x;
}
template<typename Dimension>
inline
SlideSurfaceMethod
SolidFSISPHHydroBase<Dimension>::
slideSurfaceMethod() const {
  return mSlideSurfaceMethod;
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
color() const {
  return mColor;
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
// Ref to the slide surface obj
//------------------------------------------------------------------------------
template<typename Dimension>
inline
SlideSurface<Dimension>&
SolidFSISPHHydroBase<Dimension>::slideSurface() const {
  return mSlideSurface;
}

}
