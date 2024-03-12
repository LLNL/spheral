//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScalev2
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Mon Mar 11 10:36:21 PDT 2024
//----------------------------------------------------------------------------//
#include "ASPHSmoothingScalev2.hh"
#include "Geometry/EigenStruct.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/rotationMatrix.hh"

#include <cmath>

namespace Spheral {

using std::min;
using std::max;
using std::abs;
using std::pow;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScalev2<Dimension>::
ASPHSmoothingScalev2():
  ASPHSmoothingScale<Dimension>() {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScalev2<Dimension>::
ASPHSmoothingScalev2(const ASPHSmoothingScalev2<Dimension>& rhs):
  ASPHSmoothingScale<Dimension>(rhs) {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScalev2<Dimension>&
ASPHSmoothingScalev2<Dimension>::
operator=(const ASPHSmoothingScalev2& rhs) {
  ASPHSmoothingScale<Dimension>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScalev2<Dimension>::
~ASPHSmoothingScalev2<Dimension>() {
}

//------------------------------------------------------------------------------
// Compute an idealized new H based on the given moments.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
ASPHSmoothingScalev2<Dimension>::
idealSmoothingScale(const SymTensor& H,
                    const Vector& pos,
                    const Scalar zerothMoment,
                    const SymTensor& secondMoment,
                    const TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh,
                    const ConnectivityMap<Dimension>& connectivityMap,
                    const unsigned nodeListi,
                    const unsigned i) const {

  // Pre-conditions.
  REQUIRE(H.Determinant() > 0.0);
  REQUIRE(zerothMoment >= 0.0);
  REQUIRE(secondMoment.Determinant() >= 0.0);

  // const double tiny = 1.0e-50;
  // const double tolerance = 1.0e-5;

  // If there is no information to be had (no neighbors), just double the current H vote
  // and bail
  if (secondMoment.Determinant() == 0.0) return 0.5*H;

  // Decompose the second moment tensor into it's eigen values/vectors.
  const auto Psi_eigen = secondMoment.eigenVectors();

  // Iterate over the eigen values and build the new H tensor in the kernel frame.
  SymTensor HnewInv;
  for (auto nu = 0u; nu < Dimension::nDim; ++nu) {
    const auto lambdaPsi = Psi_eigen.eigenValues(nu);
    const auto evec = Psi_eigen.eigenVectors.getColumn(nu);
    const auto h0 = 1.0/(H*evec).magnitude();

    // Query the kernel for the equivalent nodes per smoothing scale in this direction
    auto currentNodesPerSmoothingScale = W.equivalentNodesPerSmoothingScaleASPH(lambdaPsi);
    CHECK2(currentNodesPerSmoothingScale > 0.0, "Bad estimate for nPerh effective from kernel: " << currentNodesPerSmoothingScale);

    // The (limited) ratio of the desired to current nodes per smoothing scale.
    const Scalar s = min(4.0, max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30)));
    CHECK(s > 0.0);

    HnewInv(nu, nu) = h0*s;
  }

  // Rotate to the lab frame.
  const auto evec0 = Psi_eigen.eigenVectors.getColumn(0);
  const auto T = rotationMatrix(evec0).Transpose();
  HnewInv.rotationalTransform(T);

  // That's it
  return HnewInv.Inverse();
}

}
