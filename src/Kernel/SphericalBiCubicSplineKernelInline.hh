namespace Spheral {

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
SphericalBiCubicSplineKernel::operator==(const SphericalBiCubicSplineKernel& rhs) const {
  return true;
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
const TableKernel<Dim<3>>&
SphericalBiCubicSplineKernel::baseKernel3d() const {
  return mBaseKernel3d;
}

inline
const TableKernel<Dim<1>>&
SphericalBiCubicSplineKernel::baseKernel1d() const {
  return mBaseKernel1d;
}

inline
double
SphericalBiCubicSplineKernel::etamax() const {
  return metamax;
}

}
