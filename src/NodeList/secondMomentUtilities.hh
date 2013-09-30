//------------------------------------------------------------------------------
// SecondMomentUtilities
//
// A collection of standalone methods used in measure the second moment of the
// local node distribution.
//------------------------------------------------------------------------------
#ifndef __Spheral__secondMomentUtilities__
#define __Spheral__secondMomentUtilities__

//------------------------------------------------------------------------------
// Compute the symmetrized version of a tensor via geometric mean.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::SymTensor
geometricSymmetrize(const typename Dimension::Tensor& A) {

  // Compute A^2.
  typedef typename Dimension::SymTensor SymTensor;
  const SymTensor A2 = (A*A).Symmetric();

  // Now take the square root.
  const typename SymTensor::EigenStructType eigen = A2.eigenVectors();
  SymTensor result;
  for (int i = 0; i != Dimension::nDim; ++i) result(i,i) = sqrt(eigen.eigenValues(i));
  result.rotationalTransform(eigen.eigenVectors);

  return result;
}

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
inline
double
reduceInnerWeight(const double eta, const double etac) {
  REQUIRE(eta >= 0.0);
  REQUIRE(etac > 0.0);
  if (eta > etac) {
    return 1.0;
  } else {
    const double thpt = 0.5*(1.0 + sin((eta/etac - 0.5)*M_PI));
    return thpt*thpt;
  }
}

//------------------------------------------------------------------------------
// Helper to compute the weighted neighbor sum contribution for determining H.
//------------------------------------------------------------------------------
template<typename KernelType>
inline
double
computeNeighborWeight(const double& eta,
                      const KernelType& W) {
  REQUIRE(eta >= 0.0);
  const double Wi = W(eta, 1.0)*eta/(eta + 1.0e-10);
  return Wi;
}

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// void
// incrementWeightedGeometricalMean(typename Dimension::SymTensor& mean,
//                                  const typename Dimension::SymTensor& H,
//                                  const typename Dimension::Scalar& weight) {

//   // Pre-conditions.
//   REQUIRE(H.Determinant() > 0.0);
//   REQUIRE(weight >= 0.0);

//   // Determine the weighted H inverse.
//   const double nDimInv = 1.0/Dimension::nDim;
//   const typename Dimension::SymTensor Hi = H / Dimension::rootnu(H.Determinant());
//   CHECK(fuzzyEqual(Hi.Determinant(), 1.0));
//   const typename Dimension::SymTensor::EigenStructType eigen = Hi.eigenVectors();
//   typename Dimension::SymTensor wHi;
//   for (int i = 0; i != Dimension::nDim; ++i) wHi(i,i) = pow(eigen.eigenValues(i), 0.5*weight);
//   wHi.rotationalTransform(eigen.eigenVectors);

//   // Apply the weighted Hi to the cumulative result.
//   typename Dimension::SymTensor newmean = (wHi*mean*wHi).Symmetric();
//   mean = newmean;
// }

#endif
