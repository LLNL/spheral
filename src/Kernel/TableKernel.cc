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

namespace {  // anonymous

//------------------------------------------------------------------------------
// Sum the Kernel values for the given stepsize (SPH)
//------------------------------------------------------------------------------
inline
double
sumKernelValues(const TableKernel<Dim<1>>& W,
                const double deta) {
  REQUIRE(deta > 0);
  double result = 0.0;
  double etar = deta;
  while (etar < W.kernelExtent()) {
    result += 2.0*std::abs(W.gradValue(etar, 1.0));
    etar += deta;
  }
  return result;
}

inline
double
sumKernelValues(const TableKernel<Dim<2>>& W,
                const double deta) {
  REQUIRE(deta > 0);
  double result = 0.0;
  double etar = deta;
  while (etar < W.kernelExtent()) {
    result += 2.0*M_PI*etar/deta*std::abs(W.gradValue(etar, 1.0));
    etar += deta;
  }
  return sqrt(result);
}

inline
double
sumKernelValues(const TableKernel<Dim<3>>& W,
                const double deta) {
  REQUIRE(deta > 0);
  double result = 0.0;
  double etar = deta;
  while (etar < W.kernelExtent()) {
    result += 4.0*M_PI*FastMath::square(etar/deta)*std::abs(W.gradValue(etar, 1.0));
    etar += deta;
  }
  return pow(result, 1.0/3.0);
}

//------------------------------------------------------------------------------
// Sum the Kernel values for the given stepsize (ASPH)
// We do these on a lattice pattern since the coordinates of the points are
// used.
//------------------------------------------------------------------------------
inline
double
sumKernelValuesASPH(const TableKernel<Dim<1>>& W,
                    const double deta) {
  REQUIRE(deta > 0);
  auto result = 0.0;
  auto etax = deta;
  while (etax < W.kernelExtent()) {
    result += 2.0*std::abs(W.gradValue(etax, 1.0)) * etax*etax;
    etax += deta;
  }
  return result;
}

inline
double
sumKernelValuesASPH(const TableKernel<Dim<2>>& W,
                    const double deta) {
  REQUIRE(deta > 0);
  Dim<2>::SymTensor result;
  double etay = 0.0;
  while (etay < W.kernelExtent()) {
    double etax = 0.0;
    while (etax < W.kernelExtent()) {
      const Dim<2>::Vector eta(etax, etay);
      auto dresult = std::abs(W.gradValue(eta.magnitude(), 1.0)) * eta.selfdyad();
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
                    const double deta) {
  REQUIRE(deta > 0);
  Dim<3>::SymTensor result;
  double etaz = 0.0;
  while (etaz < W.kernelExtent()) {
    double etay = 0.0;
    while (etay < W.kernelExtent()) {
      double etax = 0.0;
      while (etax < W.kernelExtent()) {
        const Dim<3>::Vector eta(etax, etay, etaz);
        auto dresult = std::abs(W.gradValue(eta.magnitude(), 1.0)) * eta.selfdyad();
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

//------------------------------------------------------------------------------
// Compute the (f1,f2) integrals relation for the given zeta = r/h 
// (RZ corrections).
//------------------------------------------------------------------------------
template<typename KernelType>
class volfunc {
  const KernelType& W;
public:
  volfunc(const KernelType& W): W(W) {}
  double operator()(const double eta) const {
    return W.kernelValue(eta, 1.0);
  }
};

template<typename KernelType>
class f1func {
  const KernelType& W;
  double zeta;
public:
  f1func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    return abs(safeInvVar(zeta)*eta)*W.kernelValue(abs(zeta - eta), 1.0);
  }
};


template<typename KernelType>
class f2func {
  const KernelType& W;
  double zeta;
public:
  f2func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    return safeInvVar(zeta*zeta)*eta*abs(eta)*W.kernelValue(abs(zeta - eta), 1.0);
  }
};


template<typename KernelType>
class gradf1func {
  const KernelType& W;
  double zeta;
public:
  gradf1func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    const double Wu = W.kernelValue(abs(zeta - eta), 1.0);
    const double gWu = W.gradValue(abs(zeta - eta), 1.0);
    const double gf1inv = safeInvVar(zeta)*abs(eta)*gWu - safeInvVar(zeta*zeta)*abs(eta)*Wu;
    if (eta < 0.0) {
      return -gf1inv;
    } else {
      return gf1inv;
    }
  }
};


template<typename KernelType>
double
f1Integral(const KernelType& W,
           const double zeta,
           const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return safeInvVar(simpsonsIntegration<f1func<KernelType>, double, double>(f1func<KernelType>(W, zeta), 
                                                                            zeta - etaMax, 
                                                                            zeta + etaMax,
                                                                            numbins));
}

template<typename KernelType>
double
f2Integral(const KernelType& W,
           const double zeta,
           const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return safeInvVar(simpsonsIntegration<f2func<KernelType>, double, double>(f2func<KernelType>(W, zeta), 
                                                                            zeta - etaMax, 
                                                                            zeta + etaMax,
                                                                            numbins));
}

template<typename KernelType>
double
gradf1Integral(const KernelType& W,
               const double zeta,
               const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return simpsonsIntegration<gradf1func<KernelType>, double, double>(gradf1func<KernelType>(W, zeta), 
                                                                     zeta - etaMax, 
                                                                     zeta + etaMax,
                                                                     numbins);
}

}  // anonymous

//------------------------------------------------------------------------------
// Construct from a kernel.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename KernelType>
TableKernel<Dimension>::TableKernel(const KernelType& kernel,
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
  mWsumLookup(),
  mNperhLookupASPH(),
  mWsumLookupASPH() {

  // Gotta have a minimally reasonable nperh range
  if (mMaxNperh <= mMinNperh) mMaxNperh = 4.0*mMinNperh;

  // Pre-conditions.
  VERIFY(mNumPoints > 0);
  VERIFY(mMinNperh > 0.0 and mMaxNperh > mMinNperh);

  // Set the volume normalization and kernel extent.
  this->setVolumeNormalization(1.0); // (kernel.volumeNormalization() / Dimension::pownu(hmult));  // We now build this into the tabular kernel values.
  this->setKernelExtent(kernel.kernelExtent());
  this->setInflectionPoint(kernel.inflectionPoint());

  // Set the interpolation methods for looking up nperh
  mWsumLookup.initialize(mMinNperh, mMaxNperh, numPoints,
                         [&](const double x) -> double { return sumKernelValues(*this, 1.0/x); });
  mNperhLookup.initialize(mWsumLookup(mMinNperh), mWsumLookup(mMaxNperh), numPoints,
                          [&](const double Wsum) -> double { return bisectRoot([&](const double nperh) { return mWsumLookup(nperh) - Wsum; }, mMinNperh, mMaxNperh); });

  // ASPH variants
  mWsumLookupASPH.initialize(mMinNperh, mMaxNperh, numPoints,
                             [&](const double x) -> double { return sumKernelValuesASPH(*this, 1.0/x); });
  mNperhLookupASPH.initialize(mWsumLookupASPH(mMinNperh), mWsumLookupASPH(mMaxNperh), numPoints,
                              [&](const double Wsum) -> double { return bisectRoot([&](const double nperh) { return mWsumLookupASPH(nperh) - Wsum; }, mMinNperh, mMaxNperh); });

  // // Make nperh lookups monotonic
  // mWsumLookup.makeMonotonic();
  // mNperhLookup.makeMonotonic();
  // mWsumLookupASPH.makeMonotonic();
  // mNperhLookupASPH.makeMonotonic();
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
  mWsumLookup(rhs.mWsumLookup),
  mNperhLookupASPH(rhs.mNperhLookupASPH),
  mWsumLookupASPH(rhs.mWsumLookupASPH) {
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
    mNperhLookupASPH = rhs.mNperhLookupASPH;
    mWsumLookupASPH = rhs.mWsumLookupASPH;
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
          (mGrad2Interp == rhs.mGrad2Interp));
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
// Same as above for ASPH second moment measurement
// lambda_psi here is the 1/sqrt(eigenvalue) of the second moment tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
TableKernel<Dimension>::
equivalentNodesPerSmoothingScaleASPH(const Scalar lambdaPsi) const {
  return std::max(0.0, mNperhLookupASPH(lambdaPsi));
}

template<typename Dimension>
typename Dimension::Scalar
TableKernel<Dimension>::
equivalentLambdaPsiASPH(const Scalar nPerh) const {
  return std::max(0.0, mWsumLookupASPH(nPerh));
}

}
