namespace Spheral {

//------------------------------------------------------------------------------
// set/get the generalized density exponent (Monaghan 1992)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
alpha(double x) {
  mAlpha = x;
}


template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
alpha() const {
  return mAlpha;
}

template<typename Dimension>
inline
void
SolidFSISPHHydroBase<Dimension>::
diffusionCoefficient(double x) {
  mDiffusionCoefficient = x;
}

template<typename Dimension>
inline
double
SolidFSISPHHydroBase<Dimension>::
diffusionCoefficient() const {
  return mDiffusionCoefficient;
}



//------------------------------------------------------------------------------
// swtich to turn on density sum for different nodeLists
//------------------------------------------------------------------------------
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

}
