namespace Spheral {

//------------------------------------------------------------------------------
// set/get the generalized density exponent (Monaghan 1992)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FSISolidSPHHydroBase<Dimension>::
alpha(double x) {
  mAlpha = x;
}


template<typename Dimension>
inline
double
FSISolidSPHHydroBase<Dimension>::
alpha() const {
  return mAlpha;
}


//------------------------------------------------------------------------------
// swtich to turn on density sum for different nodeLists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FSISolidSPHHydroBase<Dimension>::
sumDensityNodeListSwitch(std::vector<int> x) {
  mSumDensityNodeListSwitch = x;
}


template<typename Dimension>
inline
std::vector<int>
FSISolidSPHHydroBase<Dimension>::
sumDensityNodeListSwitch() const {
  return mSumDensityNodeListSwitch;
}

}
