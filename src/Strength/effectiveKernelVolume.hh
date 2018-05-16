//------------------------------------------------------------------------------
// Compute the (dimension dependent) effective node volume for a given H.
//------------------------------------------------------------------------------
inline
double
effectiveKernelVolume(const Dim<1>::SymTensor& H, const double kernelExtent) {
  REQUIRE(distinctlyGreaterThan(H.Determinant(), 0.0));
  return 2.0*kernelExtent/(H.Determinant());
}

inline
double
effectiveKernelVolume(const Dim<2>::SymTensor& H, const double kernelExtent) {
  REQUIRE(distinctlyGreaterThan(H.Determinant(), 0.0));
  return M_PI*kernelExtent*kernelExtent/(H.Determinant());
}

inline
double
effectiveKernelVolume(const Dim<3>::SymTensor& H, const double kernelExtent) {
  REQUIRE(distinctlyGreaterThan(H.Determinant(), 0.0));
  return 4.0/3.0*M_PI*FastMath::cube(kernelExtent)/(H.Determinant());
}
