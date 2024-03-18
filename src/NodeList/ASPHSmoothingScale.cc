//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScale
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 13:50:49 PDT 2005
//----------------------------------------------------------------------------//
#include "ASPHSmoothingScale.hh"
#include "Geometry/EigenStruct.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/bisectRoot.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Mesh/Mesh.hh"

#include <cmath>

namespace Spheral {

using std::vector;
using std::min;
using std::max;
using std::abs;
using std::pow;

namespace {

//------------------------------------------------------------------------------
// Convert a given number of neighbors to the equivalent 1D "radius" in nodes.
//------------------------------------------------------------------------------
template<typename Dimension> double equivalentRadius(const double n);

// 1D
template<>
inline
double
equivalentRadius<Dim<1> >(const double n) {
  return 0.5*n;
}

// 2D
template<>
inline
double
equivalentRadius<Dim<2> >(const double n) {
  return std::sqrt(n/M_PI);
}

// 3D
template<>
inline
double
equivalentRadius<Dim<3> >(const double n) {
  return Dim<3>::rootnu(3.0*n/(4.0*M_PI));
}

//------------------------------------------------------------------------------
// Apply a distortion to a symmetric tensor, returning a new symmetric tensor.
//------------------------------------------------------------------------------
template<typename T> 
inline
T
applyStretchToHinv(const T& stretch, const T& H0) {
  const T stretch12 = stretch.sqrt();
  return (stretch12*H0*stretch12).Symmetric();
}

template<> 
inline
Dim<2>::SymTensor
applyStretchToHinv(const Dim<2>::SymTensor& stretch, const Dim<2>::SymTensor& H0) {
  typedef Dim<2>::Tensor Tensor;
  typedef Dim<2>::SymTensor SymTensor;

  const Tensor A = (H0*stretch).Inverse();
  const double Adet = A.Determinant();
  CHECK(distinctlyGreaterThan(Adet, 0.0));

  const double Bxx = A.xx()*A.xx() + A.xy()*A.xy();
  const double Bxy = A.xx()*A.yx() + A.xy()*A.yy();
  const double Byy = A.yx()*A.yx() + A.yy()*A.yy();
  const SymTensor B(Bxx, Bxy,
                    Bxy, Byy);

  SymTensor result = B.sqrt();
  const double det = result.Determinant();
  CHECK(distinctlyGreaterThan(det, 0.0));

  result *= sqrt(Adet/det);
  ENSURE(fuzzyEqual(result.Determinant(), Adet, 1.0e-10));

  return result.Inverse();
}

//------------------------------------------------------------------------------
// Compute a weight representing how different two H tensors are.
// result = 1 if they're the same, \in [0, 1[ if they're different.
//------------------------------------------------------------------------------
template<typename Dimension> 
inline
double
Hdifference(const typename Dimension::SymTensor& H1,
            const typename Dimension::SymTensor& H2) {

  // Pre-conditions.
  REQUIRE(fuzzyEqual(H1.Determinant(), 1.0, 1.0e-5));
  REQUIRE(fuzzyEqual(H2.Determinant(), 1.0, 1.0e-5));

  typedef typename Dimension::Tensor Tensor;

  const Tensor K = H1.Inverse() * H2 - Tensor::one;
  const double result = 1.0 - min(1.0, K.doubledot(K));

  ENSURE(result >= 0.0 && result <= 1.0);
  return result;
}

//------------------------------------------------------------------------------
// Compute the new symmetric H inverse from the non-symmetric "A" tensor.
//------------------------------------------------------------------------------
inline
Dim<1>::SymTensor
computeHinvFromA(const Dim<1>::Tensor&) {
  return Dim<1>::SymTensor::one;
}

inline
Dim<2>::SymTensor
computeHinvFromA(const Dim<2>::Tensor& A) {
  REQUIRE(fuzzyEqual(A.Determinant(), 1.0, 1.0e-8));
  typedef Dim<2>::SymTensor SymTensor;
  const double A11 = A.xx();
  const double A12 = A.xy();
  const double A21 = A.yx();
  const double A22 = A.yy();
  const SymTensor Hinv2(A11*A11 + A12*A12, A11*A21 + A12*A22,
                        A11*A21 + A12*A22, A21*A21 + A22*A22);
  CHECK(distinctlyGreaterThan(Hinv2.Determinant(), 0.0));
  const SymTensor result = Hinv2.sqrt() / sqrt(sqrt(Hinv2.Determinant()));
  ENSURE(fuzzyEqual(result.Determinant(), 1.0, 1.0e-8));
  return result;
}

inline
Dim<3>::SymTensor
computeHinvFromA(const Dim<3>::Tensor&) {
  return Dim<3>::SymTensor::one;
}

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

}  // anonymous namespace

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScale<Dimension>::
ASPHSmoothingScale(const TableKernel<Dimension>& W,
                   const Scalar targetNperh,
                   const size_t numPoints):
  SmoothingScaleBase<Dimension>(),
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
ASPHSmoothingScale<Dimension>::
ASPHSmoothingScale(const ASPHSmoothingScale<Dimension>& rhs):
  SmoothingScaleBase<Dimension>(rhs),
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
ASPHSmoothingScale<Dimension>&
ASPHSmoothingScale<Dimension>::
operator=(const ASPHSmoothingScale& rhs) {
  SmoothingScaleBase<Dimension>::operator=(rhs);
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
ASPHSmoothingScale<Dimension>::
~ASPHSmoothingScale<Dimension>() {
}

//------------------------------------------------------------------------------
// Time derivative of the smoothing scale.
//------------------------------------------------------------------------------
// 1-D case same as SPH.
#ifdef SPHERAL1DINSTANTIATION
template<>
Dim<1>::SymTensor
ASPHSmoothingScale<Dim<1> >::
smoothingScaleDerivative(const Dim<1>::SymTensor& H,
                         const Dim<1>::Vector& /*pos*/,
                         const Dim<1>::Tensor& DvDx,
                         const Dim<1>::Scalar /*hmin*/,
                         const Dim<1>::Scalar /*hmax*/,
                         const Dim<1>::Scalar /*hminratio*/,
                         const Dim<1>::Scalar /*nPerh*/) const {
  return -H*DvDx.Trace();
}
#endif

#ifdef SPHERAL2DINSTANTIATION
// 2-D ASPH tensor evolution.
template<>
Dim<2>::SymTensor
ASPHSmoothingScale<Dim<2> >::
smoothingScaleDerivative(const Dim<2>::SymTensor& H,
                         const Dim<2>::Vector& /*pos*/,
                         const Dim<2>::Tensor& DvDx,
                         const Dim<2>::Scalar /*hmin*/,
                         const Dim<2>::Scalar /*hmax*/,
                         const Dim<2>::Scalar /*hminratio*/,
                         const Dim<2>::Scalar /*nPerh*/) const {
  REQUIRE(H.Trace() > 0.0);
  const Scalar thetaDot = 
    (H.xx()*DvDx.xy() - H.yy()*DvDx.yx() - H.yx()*(DvDx.xx() - DvDx.yy()))/
    H.Trace();
  SymTensor result;
  result.xx(H.yx()*(thetaDot - DvDx.yx()) - H.xx()*DvDx.xx());
  result.yx(-(H.xx()*thetaDot + H.yx()*DvDx.xx() + H.yy()*DvDx.yx()));
  result.yy(-H.yx()*(thetaDot + DvDx.xy()) - H.yy()*DvDx.yy());
  return result;
}
#endif

#ifdef SPHERAL3DINSTANTIATION
// 3-D ASPH tensor evolution.
template<>
Dim<3>::SymTensor
ASPHSmoothingScale<Dim<3> >::
smoothingScaleDerivative(const Dim<3>::SymTensor& H,
                         const Dim<3>::Vector& /*pos*/,
                         const Dim<3>::Tensor& DvDx,
                         const Dim<3>::Scalar /*hmin*/,
                         const Dim<3>::Scalar /*hmax*/,
                         const Dim<3>::Scalar /*hminratio*/,
                         const Dim<3>::Scalar /*nPerh*/) const {
  REQUIRE(H.Trace() > 0.0);
  const double AA = H.xx()*DvDx.xy() - H.xy()*(DvDx.xx() - DvDx.yy()) + 
    H.xz()*DvDx.zy() - H.yy()*DvDx.yx() - H.yz()*DvDx.zx();
  const double BB = H.xx()*DvDx.xz() + H.xy()*DvDx.yz() - 
    H.xz()*(DvDx.xx() - DvDx.zz()) - H.yz()*DvDx.yx() - H.zz()*DvDx.zx();
  const double CC = H.xy()*DvDx.xz() + H.yy()*DvDx.yz() - 
    H.yz()*(DvDx.yy() - DvDx.zz()) - H.xz()*DvDx.xy() - H.zz()*DvDx.zy();
  const double thpt = H.yy() + H.zz();
  const double Ga = (H.xx() + H.yy())*thpt - H.xz()*H.xz();
  const double Gb = (H.yy() + H.zz())*H.yz() + H.xy()*H.xz();
  const double Gc = (H.xx() + H.zz())*thpt - H.xy()*H.xy();
  const double Gd = thpt*AA + H.xz()*CC;
  const double Ge = thpt*BB - H.xy()*CC;
  const double ack = 1.0/(Ga*Gc - Gb*Gb);
  const double Gdot = (Gc*Gd - Gb*Ge)*ack;
  const double Tdot = (Gb*Gd - Ga*Ge)*ack;
  const double Phidot = (H.xz()*Gdot + H.xy()*Tdot + CC)/thpt;
  SymTensor result;
  result.xx(-H.xx()*DvDx.xx() + H.xy()*(Gdot - DvDx.yx()) - H.xz()*(Tdot + DvDx.zx()));
  result.xy(H.yy()*Gdot - H.yz()*Tdot - H.xx()*DvDx.xy() - H.xy()*DvDx.yy() - H.xz()*DvDx.zy());
  result.xz(H.yz()*Gdot - H.zz()*Tdot - H.xx()*DvDx.xz() - H.xy()*DvDx.yz() - H.xz()*DvDx.zz());
  result.yy(H.yz()*(Phidot - DvDx.zy()) - H.xy()*(Gdot + DvDx.xy()) - H.yy()*DvDx.yy());
  result.yz(H.xy()*Tdot - H.yy()*Phidot - H.xz()*DvDx.xy() - H.yz()*DvDx.yy() - H.zz()*DvDx.zy());
  result.zz(H.xz()*(Tdot - DvDx.xz()) - H.yz()*(Phidot + DvDx.yz()) - H.zz()*DvDx.zz());
  return result;
}
#endif
  
//------------------------------------------------------------------------------
// Compute an idealized new H based on the given moments.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
ASPHSmoothingScale<Dimension>::
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
// Determine a new smoothing scale as a replacement for the old, using assorted
// limiting on the ideal H measurement.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
ASPHSmoothingScale<Dimension>::
newSmoothingScale(const SymTensor& H,
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

  const double tolerance = 1.0e-5;

  // Get the ideal H vote.
  const SymTensor Hideal = idealSmoothingScale(H, 
                                               pos,
                                               zerothMoment,
                                               firstMoment,
                                               secondMomentEta,
                                               secondMomentLab,
                                               W,
                                               hmin,
                                               hmax,
                                               hminratio,
                                               nPerh,
                                               connectivityMap,
                                               nodeListi,
                                               i);

  const double Hidealscale = Dimension::rootnu(Hideal.Determinant());
  const SymTensor Hidealhatinv = Hideal.Inverse() * Hidealscale;
  CONTRACT_VAR(tolerance);
  CHECK(fuzzyEqual(Hidealhatinv.Determinant(), 1.0, tolerance));

  // Compute a weighting factor measuring how different the target H is from the old.
  const SymTensor H0hatinv = H.Inverse() * Dimension::rootnu(H.Determinant());
  CHECK(fuzzyEqual(H0hatinv.Determinant(), 1.0, tolerance));
  const Scalar st = sqrt(Hdifference<Dimension>(Hidealhatinv, H0hatinv));
  CHECK(st >= 0.0 && st <= 1.0);

  // Geometrically combine the old shape with the ideal.
  const double w1 = 0.4*(1.0 + st);
  const double w0 = 1.0 - w1;
  CHECK(w0 >= 0.0 && w0 <= 1.0);
  CHECK(w1 >= 0.0 && w1 <= 1.0);
  CHECK(fuzzyEqual(w0 + w1, 1.0));
  const typename SymTensor::EigenStructType eigen0 = H0hatinv.eigenVectors();
  const typename SymTensor::EigenStructType eigen1 = Hidealhatinv.eigenVectors();
  CHECK(eigen0.eigenValues.minElement() > 0.0);
  CHECK(eigen1.eigenValues.minElement() > 0.0);
  SymTensor wH0 = constructSymTensorWithPowDiagonal(eigen0.eigenValues, w0);
  SymTensor wH1 = constructSymTensorWithPowDiagonal(eigen1.eigenValues, 0.5*w1);
  wH0.rotationalTransform(eigen0.eigenVectors);
  wH1.rotationalTransform(eigen1.eigenVectors);
  SymTensor H1hatinv = (wH1*wH0*wH1).Symmetric();
  CHECK(H1hatinv.Determinant() > 0.0);
  H1hatinv /= Dimension::rootnu(H1hatinv.Determinant());
  CONTRACT_VAR(tolerance);
  CHECK(fuzzyEqual(H1hatinv.Determinant(), 1.0, tolerance));

  // Scale the answer to recover the determinant.
  const SymTensor H1inv = H1hatinv/Hidealscale;

  // Apply limiting to build our final answer.
  const typename SymTensor::EigenStructType eigen = H1inv.eigenVectors();
  const double effectivehmin = max(hmin,
                                   hminratio*min(hmax, eigen.eigenValues.maxElement()));
  CHECK(effectivehmin >= hmin && effectivehmin <= hmax);
  CHECK(fuzzyGreaterThanOrEqual(effectivehmin/min(hmax, eigen.eigenValues.maxElement()), hminratio));
  SymTensor result;
  for (int i = 0; i != Dimension::nDim; ++i) result(i,i) = 1.0/max(effectivehmin, min(hmax, eigen.eigenValues(i)));
  result.rotationalTransform(eigen.eigenVectors);

  // We're done!
  BEGIN_CONTRACT_SCOPE
  {
    const Vector eigenValues = result.eigenValues();
    ENSURE(distinctlyGreaterThan(eigenValues.minElement(), 0.0));
    ENSURE(fuzzyGreaterThanOrEqual(1.0/eigenValues.maxElement(), hmin, 1.0e-5));
    ENSURE(fuzzyLessThanOrEqual(1.0/eigenValues.minElement(), hmax, 1.0e-5));
    ENSURE2(fuzzyGreaterThanOrEqual(eigenValues.minElement()/eigenValues.maxElement(), hminratio, 1.e-3), (eigenValues.minElement()/eigenValues.maxElement()) << " " << hminratio);
  }
  END_CONTRACT_SCOPE

  return result;

}

//------------------------------------------------------------------------------
// Use the volumes of tessellation to set the new Hs.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
ASPHSmoothingScale<Dimension>::
idealSmoothingScale(const SymTensor& /*H*/,
                    const Mesh<Dimension>& mesh,
                    const typename Mesh<Dimension>::Zone& zone,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh) const {

  const vector<unsigned>& nodeIDs = zone.nodeIDs();
  const Scalar vol = zone.volume();
  const Vector zc = zone.position();
  CONTRACT_VAR(vol);
  CHECK(vol > 0.0);

  // Measure the second moment of the zone shape.
  SymTensor psi;
  for (unsigned j = 0; j != nodeIDs.size(); ++j) {
    const Vector dn = mesh.node(nodeIDs[j]).position() - zc;
    psi += dn.selfdyad();
  }
  psi /= nodeIDs.size();

  // Take the square root to get the shape, then rescale to get
  // the correct volume.
  SymTensor H1inv = psi.sqrt();
  H1inv *= 2.0*nPerh;

  // Apply limits.
  const typename SymTensor::EigenStructType eigen = H1inv.eigenVectors();
  const double effectivehmin = max(hmin,
                                   hminratio*min(hmax, eigen.eigenValues.maxElement()));
  CHECK(effectivehmin >= hmin && effectivehmin <= hmax);
  CHECK(fuzzyGreaterThanOrEqual(effectivehmin/min(hmax, eigen.eigenValues.maxElement()), hminratio));
  SymTensor result;
  for (unsigned j = 0; j != Dimension::nDim; ++j) {
    result(j,j) = 1.0/max(effectivehmin, min(hmax, eigen.eigenValues(j)));
  }
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// Determine the number of nodes per smoothing scale implied by the given
// sum of kernel values.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ASPHSmoothingScale<Dimension>::
equivalentNodesPerSmoothingScale(const Scalar lambdaPsi) const {
  return std::max(0.0, mNperhLookup(lambdaPsi));
}

//------------------------------------------------------------------------------
// Determine the effective Wsum we would expect for the given n per h.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ASPHSmoothingScale<Dimension>::
equivalentLambdaPsi(const Scalar nPerh) const {
  return std::max(0.0, mWsumLookup(nPerh));
}

}
