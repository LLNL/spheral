//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#include "Eigen/Dense"

#include "TableKernel.hh"

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/bisectRoot.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/safeInv.hh"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct from a kernel.
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>::TableKernel(const TableKernel<Dimension>& kernel,
                                    const unsigned numPoints,
                                    const typename Dimension::Scalar minNperh,
                                    const typename Dimension::Scalar maxNperh):
  Kernel<Dimension, TableKernel<Dimension>>(),
  mNumPoints(numPoints),
  mMinNperh(std::max(minNperh, 1.1/kernel.kernelExtent())),
  mMaxNperh(maxNperh),
  mInterp(0.0, kernel.kernelExtent(), numPoints,      [&](const double x) { return kernel(x, 1.0); }),
  mGradInterp(0.0, kernel.kernelExtent(), numPoints,  [&](const double x) { return kernel.grad(x, 1.0); }),
  mGrad2Interp(0.0, kernel.kernelExtent(), numPoints, [&](const double x) { return kernel.grad2(x, 1.0); }),
  mNperhLookup(),
  mWsumLookup() {

  // Gotta have a minimally reasonable nperh range
  if (mMaxNperh <= mMinNperh) mMaxNperh = 4.0*mMinNperh;

  // Pre-conditions.
  VERIFY2(mNumPoints > 0, "TableKernel ERROR: require numPoints > 0 : " << mNumPoints);
  VERIFY2(mMinNperh > 0.0 and mMaxNperh > mMinNperh, "TableKernel ERROR: Bad (minNperh, maxNperh) range: (" << mMinNperh << ", " << mMaxNperh << ")");

  // Set the volume normalization and kernel extent.
  this->setVolumeNormalization(1.0); // (kernel.volumeNormalization() / Dimension::pownu(hmult));  // We now build this into the tabular kernel values.
  this->setKernelExtent(kernel.kernelExtent());
  this->setInflectionPoint(kernel.inflectionPoint());

  // Set the interpolation methods for looking up nperh (SPH methodology)
  mWsumLookup.initialize(mMinNperh, mMaxNperh, numPoints,
                         [&](const double x) -> double { return sumKernelValues(*this, x); });
  mNperhLookup.initialize(mWsumLookup(mMinNperh), mWsumLookup(mMaxNperh), numPoints,
                          [&](const double Wsum) -> double { return bisectRoot([&](const double nperh) { return mWsumLookup(nperh) - Wsum; }, mMinNperh, mMaxNperh); });

  // Make nperh lookups monotonic
  mWsumLookup.makeMonotonic();
  mNperhLookup.makeMonotonic();
}

//------------------------------------------------------------------------------
// Copy
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>::
TableKernel(const TableKernel<Dimension>& rhs):
  Kernel<Dimension, TableKernel<Dimension>>(rhs),
  mNumPoints(rhs.mNumPoints),
  mMinNperh(rhs.mMinNperh),
  mMaxNperh(rhs.mMaxNperh),
  mInterp(rhs.mInterp),
  mGradInterp(rhs.mGradInterp),
  mGrad2Interp(rhs.mGrad2Interp),
  mNperhLookup(rhs.mNperhLookup),
  mWsumLookup(rhs.mWsumLookup) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>::
~TableKernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>&
TableKernel<Dimension>::
operator=(const TableKernel<Dimension>& rhs) {
  if (this != &rhs) {
    Kernel<Dimension, TableKernel<Dimension>>::operator=(rhs);
    mNumPoints = rhs.mNumPoints;
    mMinNperh = rhs.mMinNperh;
    mMaxNperh = rhs.mMaxNperh;
    mInterp = rhs.mInterp;
    mGradInterp = rhs.mGradInterp;
    mGrad2Interp = rhs.mGrad2Interp;
    mNperhLookup = rhs.mNperhLookup;
    mWsumLookup = rhs.mWsumLookup;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TableKernel<Dimension>::
operator==(const TableKernel<Dimension>& rhs) const {
  return ((mInterp == rhs.mInterp) and
          (mGradInterp == rhs.mGradInterp) and
          (mGrad2Interp == rhs.mGrad2Interp) and
          (mNperhLookup == rhs.mNperhLookup) and
          (mWsumLookup == rhs.mWsumLookup));
}

//------------------------------------------------------------------------------
// Kernel value for SPH smoothing scale nperh lookups
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TableKernel<Dimension>::kernelValueSPH(const Scalar etaij) const {
  REQUIRE(etaij >= 0.0);
  if (etaij < this->mKernelExtent) {
    return std::abs(mGradInterp(etaij));
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Kernel value for ASPH smoothing scale nperh lookups
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TableKernel<Dimension>::kernelValueASPH(const Scalar etaij, const Scalar nPerh) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(nPerh > 0.0);
  if (etaij < mKernelExtent) {
    const auto deta = 2.0/std::max(2.0, nPerh);
    const auto eta0 = std::max(0.0, 0.5*(mKernelExtent - deta));
    const auto eta1 = std::min(mKernelExtent, eta0 + deta);
    return (etaij <= eta0 or etaij >= eta1 ?
            0.0 :
            kernelValueSPH((etaij - eta0)/deta));
            // FastMath::square(sin(M_PI*(etaij - eta0)/deta)));
    // return std::abs(mGradInterp(etaij * std::max(1.0, 0.5*nPerh*mKernelExtent))); // * FastMath::square(sin(nPerh*M_PI*etaij));
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Determine the number of nodes per smoothing scale implied by the given
// sum of kernel values (SPH round tensor definition).
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TableKernel<Dimension>::
equivalentNodesPerSmoothingScale(const Scalar Wsum) const {
  return std::max(0.0, mNperhLookup(Wsum));
}

//------------------------------------------------------------------------------
// Determine the effective Wsum we would expect for the given n per h.
// (SPH round tensor definition).
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TableKernel<Dimension>::
equivalentWsum(const Scalar nPerh) const {
  return std::max(0.0, mWsumLookup(nPerh));
}

//------------------------------------------------------------------------------
// Look up the kernel and first derivative for a set.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TableKernel<Dimension>::kernelAndGradValues(const Scalar* etaijs,
                                            const Scalar* Hdets,
                                            Scalar* kernelValues,
                                            Scalar* gradValues,
                                            const size_t n) const {
  for (size_t i = 0; i < n; ++i) {
    kernelAndGradValue(etaijs[i], Hdets[i], kernelValues[i], gradValues[i]);
  }
}

}