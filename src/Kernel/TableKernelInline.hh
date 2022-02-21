#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"
#include "VolumeIntegrationFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::kernelValue(const double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaij < this->mKernelExtent) {
    return Hdet*mInterp(etaij);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::gradValue(const double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaij < this->mKernelExtent) {
    return Hdet*mGradInterp(etaij);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::grad2Value(const double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaij < this->mKernelExtent) {
    return Hdet*mGrad2Interp(etaij);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the kernel and gradient for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TableKernel<Dimension>::kernelAndGrad(const typename Dimension::Vector& etaj,
                                      const typename Dimension::Vector& etai,
                                      const typename Dimension::SymTensor& H,
                                      typename Dimension::Scalar& W,
                                      typename Dimension::Vector& gradW,
                                      typename Dimension::Scalar& deltaWsum) const {
  const auto etaji = etaj - etai;
  const auto etajiMag = etaji.magnitude();
  const auto Hdet = H.Determinant();
  if (etajiMag < this->mKernelExtent) {
    const auto i0 = mInterp.lowerBound(etajiMag);
    W = Hdet*mInterp(etajiMag, i0);
    deltaWsum = Hdet*mGradInterp(etajiMag, i0);
    gradW = H*etaji.unitVector()*deltaWsum;
  } else {
    W = 0.0;
    deltaWsum = 0.0;
    gradW.Zero();
  }
}

//------------------------------------------------------------------------------
// Return the kernel and gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::pair<double, double>
TableKernel<Dimension>::kernelAndGradValue(const double etaij, const double Hdet) const {
  REQUIRE(etaij >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaij < this->mKernelExtent) {
    const auto i0 = mInterp.lowerBound(etaij);
    return std::make_pair(Hdet*mInterp(etaij, i0),
                          Hdet*mGradInterp(etaij, i0));
  } else {
    return std::make_pair(0.0, 0.0);
  }
}

//------------------------------------------------------------------------------
// Return the kernel and gradient values for a set of normalized distances.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TableKernel<Dimension>::kernelAndGradValues(const std::vector<double>& etaijs,
                                            const std::vector<double>& Hdets,
                                            std::vector<double>& kernelValues,
                                            std::vector<double>& gradValues) const {
  // Preconditions.
  const auto n = etaijs.size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(Hdets.size() == n);
    for (auto i = 0u; i < n; ++i) {
      REQUIRE(etaijs[i] >= 0.0);
      REQUIRE(Hdets[i] >= 0.0);
    }
  }
  END_CONTRACT_SCOPE

  // Prepare the results.
  kernelValues.resize(n);
  gradValues.resize(n);

  // Fill those suckers in.
  for (auto i = 0u; i < n; ++i) {
    const auto i0 = mInterp.lowerBound(etaijs[i]);
    kernelValues[i] = Hdets[i]*mInterp(etaijs[i], i0);
    gradValues[i] = Hdets[i]*mGradInterp(etaijs[i], i0);
  }
}

//------------------------------------------------------------------------------
// Return the f1 RZ correctiong for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::f1(const double /*etaij*/) const {
  VERIFY2(false, "TableKernel::f1 lookup only valid for 2D kernels.");
}

template<>
inline
double
TableKernel<Dim<2> >::f1(const double etaij) const {
  REQUIRE(etaij >= 0.0);
  if (etaij < this->mKernelExtent) {
    return mf1Interp(etaij);
  } else {
    return 1.0;
  }
}

//------------------------------------------------------------------------------
// Return the f2 RZ correctiong for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::f2(const double /*etaij*/) const {
  VERIFY2(false, "TableKernel::f2 lookup only valid for 2D kernels.");
}

template<>
inline
double
TableKernel<Dim<2> >::f2(const double etaij) const {
  REQUIRE(etaij >= 0.0);
  return mf2Interp(etaij);
}

//------------------------------------------------------------------------------
// Return the grad_r f1 RZ correctiong for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::gradf1(const double /*etaij*/) const {
  VERIFY2(false, "TableKernel::gradf1 lookup only valid for 2D kernels.");
}

template<>
inline
double
TableKernel<Dim<2> >::gradf1(const double etaij) const {
  REQUIRE(etaij >= 0.0);
  return mf1Interp.prime(etaij);
}

//------------------------------------------------------------------------------
// Return the grad_r f2 RZ correctiong for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::gradf2(const double /*etaij*/) const {
  VERIFY2(false, "TableKernel::gradf2 lookup only valid for 2D kernels.");
}

template<>
inline
double
TableKernel<Dim<2> >::gradf2(const double etaij) const {
  REQUIRE(etaij >= 0.0);
  return mf2Interp.prime(etaij);
}

//------------------------------------------------------------------------------
// Return the (f1, f2, gradf1, gradf2) RZ values for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TableKernel<Dimension>::f1Andf2(const double /*etaij*/,
                                double& /*f1*/,
                                double& /*f2*/,
                                double& /*gradf1*/,
                                double& /*gradf2*/) const {
  VERIFY2(false, "TableKernel::f1Andf2 lookup only valid for 2D kernels.");
}

template<>
inline
void
TableKernel<Dim<2> >::f1Andf2(const double etaij,
                              double& f1,
                              double& f2,
                              double& gradf1,
                              double& gradf2) const {
  REQUIRE(etaij >= 0.0);
  if (etaij < this->mKernelExtent) {
    f1 = mf1Interp(etaij);
    f2 = mf2Interp(etaij);
    gradf1 = mf1Interp.prime(etaij);
    gradf2 = mf2Interp.prime(etaij);
  } else {
    f1 = 1.0;
    f2 = 1.0;
    gradf1 = 0.0;
    gradf2 = 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the assorted tabular lookup data.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
nperhValues() const {
  return mNperhValues;
}

template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
WsumValues() const {
  return mWsumValues;
}

template<typename Dimension>
inline
size_t
TableKernel<Dimension>::
numPoints() const {
  return mNumPoints;
}

}

