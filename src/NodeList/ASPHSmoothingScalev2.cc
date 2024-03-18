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
#include "Utilities/bisectRoot.hh"

#include <cmath>

namespace Spheral {

using std::min;
using std::max;
using std::abs;
using std::pow;

namespace {

//------------------------------------------------------------------------------
// Sum the Kernel values for the given stepsize (ASPH)
// We do these on a lattice pattern since the coordinates of the points are
// used.
//------------------------------------------------------------------------------
inline
double
sumKernelValuesASPH(const TableKernel<Dim<1>>& W,
                    const double targetNperh,
                    const double nPerh) {
  REQUIRE(nPerh > 0.0);
  const auto deta = 1.0/nPerh;
  auto result = 0.0;
  auto etax = deta;
  while (etax < W.kernelExtent()) {
    result += 2.0*W.kernelValueASPH(etax, targetNperh) * etax*etax;
    etax += deta;
  }
  return result;
}

inline
double
sumKernelValuesASPH(const TableKernel<Dim<2>>& W,
                    const double targetNperh,
                    const double nPerh) {
  REQUIRE(nPerh > 0.0);
  const auto deta = 1.0/nPerh;
  Dim<2>::SymTensor result;
  double etay = 0.0;
  while (etay < W.kernelExtent()) {
    double etax = 0.0;
    while (etax < W.kernelExtent()) {
      const Dim<2>::Vector eta(etax, etay);
      auto dresult = W.kernelValueASPH(eta.magnitude(), targetNperh) * eta.selfdyad();
      if (distinctlyGreaterThan(etax, 0.0)) dresult *= 2.0;
      if (distinctlyGreaterThan(etay, 0.0)) dresult *= 2.0;
      result += dresult;
      etax += deta;
    }
    etay += deta;
  }
  const auto lambda = 0.5*(result.eigenValues().sumElements());
  return std::sqrt(lambda);
}

inline
double
sumKernelValuesASPH(const TableKernel<Dim<3>>& W,
                    const double targetNperh,
                    const double nPerh) {
  REQUIRE(nPerh > 0.0);
  const auto deta = 1.0/nPerh;
  Dim<3>::SymTensor result;
  double etaz = 0.0;
  while (etaz < W.kernelExtent()) {
    double etay = 0.0;
    while (etay < W.kernelExtent()) {
      double etax = 0.0;
      while (etax < W.kernelExtent()) {
        const Dim<3>::Vector eta(etax, etay, etaz);
        auto dresult = W.kernelValueASPH(eta.magnitude(), targetNperh) * eta.selfdyad();
        if (distinctlyGreaterThan(etax, 0.0)) dresult *= 2.0;
        if (distinctlyGreaterThan(etay, 0.0)) dresult *= 2.0;
        if (distinctlyGreaterThan(etaz, 0.0)) dresult *= 2.0;
        result += dresult;
        etax += deta;
      }
      etay += deta;
    }
    etaz += deta;
  }
  const auto lambda = (result.eigenValues().sumElements())/3.0;
  return pow(lambda, 1.0/3.0);
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScalev2<Dimension>::
ASPHSmoothingScalev2(const TableKernel<Dimension>& W,
                     const Scalar targetNperh,
                     const size_t numPoints):
  ASPHSmoothingScale<Dimension>(),
  mTargetNperh(targetNperh),
  mMinNperh(W.minNperhLookup()),
  mMaxNperh(W.maxNperhLookup()),
  mNperhLookup(),
  mWsumLookup() {

  // Preconditions
  VERIFY2(mTargetNperh >= mMinNperh, "ASPHSmoothingScale ERROR: targetNperh not in (minNperh, maxNperh) : " << mTargetNperh << " : (" << mMinNperh << ", " << mMaxNperh << ")");

  // Initalize the lookup tables for finding the effective n per h
  const auto n = numPoints > 0u ? numPoints : W.numPoints();
  mWsumLookup.initialize(mMinNperh, mMaxNperh, n,
                         [&](const double x) -> double { return sumKernelValuesASPH(W, mTargetNperh, x); });
  mNperhLookup.initialize(mWsumLookup(mMinNperh), mWsumLookup(mMaxNperh), n,
                          [&](const double Wsum) -> double { return bisectRoot([&](const double nperh) { return mWsumLookup(nperh) - Wsum; }, mMinNperh, mMaxNperh); });

  mWsumLookup.makeMonotonic();
  mNperhLookup.makeMonotonic();
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScalev2<Dimension>::
ASPHSmoothingScalev2(const ASPHSmoothingScalev2<Dimension>& rhs):
  ASPHSmoothingScale<Dimension>(rhs),
  mTargetNperh(rhs.mTargetNperh),
  mMinNperh(rhs.mMinNperh),
  mMaxNperh(rhs.mMaxNperh),
  mNperhLookup(rhs.mNperhLookup),
  mWsumLookup(rhs.mWsumLookup) {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScalev2<Dimension>&
ASPHSmoothingScalev2<Dimension>::
operator=(const ASPHSmoothingScalev2& rhs) {
  ASPHSmoothingScale<Dimension>::operator=(rhs);
  mTargetNperh = rhs.mTargetNperh;
  mMinNperh = rhs.mMinNperh;
  mMaxNperh = rhs.mMaxNperh;
  mNperhLookup = rhs.mNperhLookup;
  mWsumLookup = rhs.mWsumLookup;
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
                    const Vector& firstMoment,
                    const SymTensor& secondMomentEta,
                    const SymTensor& secondMomentLab,
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
  REQUIRE(secondMomentEta.Determinant() >= 0.0);

  // const double tiny = 1.0e-50;
  // const double tolerance = 1.0e-5;

  // If there is no information to be had (no neighbors), just double the current H vote
  // and bail
  if (secondMomentEta.Determinant() == 0.0) return 0.5*H;

  // Decompose the second moment tensor into it's eigen values/vectors.
  const auto Psi_eigen = secondMomentEta.eigenVectors();

  // Iterate over the eigen values and build the new H tensor in the kernel frame.
  SymTensor HnewInv;
  for (auto nu = 0u; nu < Dimension::nDim; ++nu) {
    const auto lambdaPsi = Psi_eigen.eigenValues(nu);
    const auto evec = Psi_eigen.eigenVectors.getColumn(nu);
    const auto h0 = 1.0/(H*evec).magnitude();

    // Query the kernel for the equivalent nodes per smoothing scale in this direction
    auto currentNodesPerSmoothingScale = this->equivalentNodesPerSmoothingScale(lambdaPsi);
    CHECK2(currentNodesPerSmoothingScale > 0.0, "Bad estimate for nPerh effective from kernel: " << currentNodesPerSmoothingScale);

    // The (limited) ratio of the desired to current nodes per smoothing scale.
    const Scalar s = min(4.0, max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30)));
    CHECK(s > 0.0);

    HnewInv(nu, nu) = h0*s;
  }

  // Rotate to the lab frame.
  HnewInv.rotationalTransform(Psi_eigen.eigenVectors);

  // That's it
  return HnewInv.Inverse();
}

//------------------------------------------------------------------------------
// Determine the number of nodes per smoothing scale implied by the given
// sum of kernel values.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ASPHSmoothingScalev2<Dimension>::
equivalentNodesPerSmoothingScale(const Scalar lambdaPsi) const {
  return std::max(0.0, mNperhLookup(lambdaPsi));
}

//------------------------------------------------------------------------------
// Determine the effective Wsum we would expect for the given n per h.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ASPHSmoothingScalev2<Dimension>::
equivalentLambdaPsi(const Scalar nPerh) const {
  return std::max(0.0, mWsumLookup(nPerh));
}

}
