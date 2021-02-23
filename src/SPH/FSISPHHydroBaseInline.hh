namespace Spheral {

//------------------------------------------------------------------------------
// set/get the generalized density exponent (Monaghan 1992)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FSISPHHydroBase<Dimension>::
alpha(double x) {
  mAlpha = x;
}


template<typename Dimension>
inline
double
FSISPHHydroBase<Dimension>::
alpha() const {
  return mAlpha;
}

template<typename Dimension>
inline
void
FSISPHHydroBase<Dimension>::
diffusionCoefficient(double x) {
  mDiffusionCoefficient = x;
}

template<typename Dimension>
inline
double
FSISPHHydroBase<Dimension>::
diffusionCoefficient() const {
  return mDiffusionCoefficient;
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
