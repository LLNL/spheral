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

template<typename Dimension>
inline
void
FSISolidSPHHydroBase<Dimension>::
diffusionCoefficient(double x) {
  mDiffusionCoefficient = x;
}

template<typename Dimension>
inline
double
FSISolidSPHHydroBase<Dimension>::
diffusionCoefficient() const {
  return mDiffusionCoefficient;
}


template<typename Dimension>
inline
void
FSISolidSPHHydroBase<Dimension>::
interfaceMethod(int x) {
  mInterfaceMethod = x;
}

template<typename Dimension>
inline
int
FSISolidSPHHydroBase<Dimension>::
interfaceMethod() const {
  return mInterfaceMethod;
}


//------------------------------------------------------------------------------
// swtich to decouple DrhoDt and DepsDt from other nodeLists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FSISolidSPHHydroBase<Dimension>::
decoupledNodeLists(std::vector<int> x) {
  mDecoupledNodeLists = x;
}


template<typename Dimension>
inline
std::vector<int>
FSISolidSPHHydroBase<Dimension>::
decoupledNodeLists() const {
  return mDecoupledNodeLists;
}

//------------------------------------------------------------------------------
// swtich to turn on density sum for different nodeLists
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FSISolidSPHHydroBase<Dimension>::
sumDensityNodeLists(std::vector<int> x) {
  mSumDensityNodeLists = x;
}


template<typename Dimension>
inline
std::vector<int>
FSISolidSPHHydroBase<Dimension>::
sumDensityNodeLists() const {
  return mSumDensityNodeLists;
}

}
