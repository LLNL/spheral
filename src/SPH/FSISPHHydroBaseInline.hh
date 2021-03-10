namespace Spheral {
//------------------------------------------------------------------------------
// set/get the generalized density exponent (Monaghan 1992)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FSISPHHydroBase<Dimension>::
surfaceForceCoefficient(double x) {
  mSurfaceForceCoefficient = x;
}


template<typename Dimension>
inline
double
FSISPHHydroBase<Dimension>::
surfaceForceCoefficient() const {
  return mSurfaceForceCoefficient;
}

template<typename Dimension>
inline
void
FSISPHHydroBase<Dimension>::
densityStabilizationCoefficient(double x) {
  mDensityStabilizationCoefficient = x;
}

template<typename Dimension>
inline
double
FSISPHHydroBase<Dimension>::
densityStabilizationCoefficient() const {
  return mDensityStabilizationCoefficient;
}

template<typename Dimension>
inline
void
FSISPHHydroBase<Dimension>::
densityDiffusionCoefficient(double x) {
  mDensityDiffusionCoefficient = x;
}

template<typename Dimension>
inline
double
FSISPHHydroBase<Dimension>::
densityDiffusionCoefficient() const {
  return mDensityDiffusionCoefficient;
}

template<typename Dimension>
inline
void
FSISPHHydroBase<Dimension>::
specificThermalEnergyDiffusionCoefficient(double x) {
  mSpecificThermalEnergyDiffusionCoefficient = x;
}

template<typename Dimension>
inline
double
FSISPHHydroBase<Dimension>::
specificThermalEnergyDiffusionCoefficient() const {
  return mSpecificThermalEnergyDiffusionCoefficient;
}
//------------------------------------------------------------------------------
// swtich to turn on density sum for different nodeLists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FSISPHHydroBase<Dimension>::
sumDensityNodeLists(std::vector<int> x) {
  mSumDensityNodeLists = x;
}


template<typename Dimension>
inline
std::vector<int>
FSISPHHydroBase<Dimension>::
sumDensityNodeLists() const {
  return mSumDensityNodeLists;
}

}
