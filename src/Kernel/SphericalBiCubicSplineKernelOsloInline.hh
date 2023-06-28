namespace Spheral {

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
SphericalBiCubicSplineKernelOslo::operator==(const SphericalBiCubicSplineKernelOslo& rhs) const {
  return true;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
const TableKernel<Dim<3>>&
SphericalBiCubicSplineKernelOslo::baseKernel3d() const {
  return mBaseKernel3d;
}

inline
const TableKernel<Dim<1>>&
SphericalBiCubicSplineKernelOslo::baseKernel1d() const {
  return mBaseKernel1d;
}

inline
double
SphericalBiCubicSplineKernelOslo::etamax() const {
  return metamax;
}

}
