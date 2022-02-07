//---------------------------------Spheral++----------------------------------//
// SphericalTableKernel
//
// Take a 3D Kernel and build a specialized 1D tabulated version appropriate
// for use with the spherical SPH algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Created by JMO, Wed Dec  2 16:41:20 PST 2020
//----------------------------------------------------------------------------//

#include "SphericalTableKernel.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

namespace {  // anonymous

//------------------------------------------------------------------------------
// Functor for doing our volume integration of the kernel  
//------------------------------------------------------------------------------
struct W3S1Func {
  const TableKernel<Dim<3>>& mW;
  double metaMax;
  W3S1Func(const TableKernel<Dim<3>>& W): mW(W), metaMax(mW.kernelExtent()) {}

  // Define a local nested functor we'll use to do the volume integral for
  // an (etaj, etai) lookup.
  struct VolFunc {
    const TableKernel<Dim<3>>& mW;
    VolFunc(const TableKernel<Dim<3>>& W): mW(W) {}
    double operator()(const double q) const {
      CHECK(q >= 0.0);
      return q * mW(q, 1.0);
    }
  };

  // Now the lookup based on (etaj, etai) -- the (rprime/h, r/h) from the paper
  double operator()(const double low, const double high) const {
    if (low >= high) return 0.0;
    return simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                        low,
                                                        high,
                                                        std::max(size_t(1000), mW.numPoints()));
  }
};

struct GradW3S1Func {
  const TableKernel<Dim<3>>& mW;
  double metaMax;
  GradW3S1Func(const TableKernel<Dim<3>>& W): mW(W), metaMax(mW.kernelExtent()) {}

  // Define a local nested functor we'll use to do the volume integral for
  // an (etaj, etai) lookup.
  struct VolFunc {
    const TableKernel<Dim<3>>& mW;
    VolFunc(const TableKernel<Dim<3>>& W): mW(W) {}
    double operator()(const double q) const {
      CHECK(q >= 0.0);
      return mW(q, 1.0) - q * mW.gradValue(q, 1.0);
    }
  };

  // Now the lookup based on (etaj, etai) -- the (rprime/h, r/h) from the paper
  double operator()(const double low, const double high) const {
    if (low >= high) return 0.0;
    return simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                        low,
                                                        high,
                                                        std::max(size_t(1000), mW.numPoints()));
  }
};

}            // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
SphericalTableKernel::SphericalTableKernel(const TableKernel<Dim<3>>& kernel):
  mInterp(0.0, kernel.kernelExtent(),
          0.0, kernel.kernelExtent(),
          std::max(size_t(200), kernel.numPoints()),
          std::max(size_t(200), kernel.numPoints()),
          W3S1Func(kernel)),
  mGradInterp(0.0, kernel.kernelExtent(),
              0.0, kernel.kernelExtent(),
              std::max(size_t(200), kernel.numPoints()),
              std::max(size_t(200), kernel.numPoints()),
              GradW3S1Func(kernel)),
  mKernel(kernel),
  metamax(kernel.kernelExtent()) {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
SphericalTableKernel::SphericalTableKernel(const SphericalTableKernel& rhs):
  mInterp(rhs.mInterp),
  mGradInterp(rhs.mGradInterp),
  mKernel(rhs.mKernel),
  metamax(rhs.metamax) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalTableKernel::~SphericalTableKernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
SphericalTableKernel&
SphericalTableKernel::operator=(const SphericalTableKernel& rhs) {
  if (this != &rhs) {
    mInterp = rhs.mInterp;
    mGradInterp = rhs.mGradInterp;
    mKernel = rhs.mKernel;
    metamax = rhs.metamax;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
bool
SphericalTableKernel::valid() const {
  return true;
}

//------------------------------------------------------------------------------
// Lookup the grad kernel for (rj/h, ri/h) = (etaj, etai)
// Using the Leibniz integral rule to differentiate this integral.
//------------------------------------------------------------------------------
Dim<1>::Vector
SphericalTableKernel::grad(const Dim<1>::Vector& etaj,
                           const Dim<1>::Vector& etai,
                           const Dim<1>::Scalar  Hdeti) const {
  REQUIRE(Hdeti >= 0.0);
  const auto ei = std::max(1e-10, etai[0]);
  const auto ej = std::max(1e-10, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) return Vector::zero;
  const auto max_bound = std::min(metamax, ei + ej);
  const auto A = (ei + ej >= metamax ?
                  0.0 :
                  max_bound*mKernel.kernelValue(max_bound, 1.0));
  const auto B = min_bound*mKernel.kernelValue(min_bound, 1.0)*(ei > ej ? 1.0 : -1.0);
  return Vector(2.0*M_PI/(ei*ej)*FastMath::pow4(Hdeti)*(A - B + 0.0*mGradInterp(min_bound, max_bound) - 1.0/ei*mInterp(min_bound, max_bound)));

  // const auto A = (ei + ej >= metamax ?
  //                 0.0 :
  //                 mKernel.kernelValue(max_bound, Hdeti) - max_bound*mKernel.gradValue(max_bound, Hdeti));
  // const auto B = min_bound*mKernel.gradValue(min_bound, Hdeti)*(ei > ej ? 1.0 : -1.0);
  // return Vector(2.0*M_PI/(ei*ej)*FastMath::cube(Hdeti)*(A + B - 1.0/ei*mInterp(min_bound, max_bound)));
  // const auto A = (ei + ej >= metamax ?
  //                 0.0 :
  //                 max_bound*mKernel.kernelValue(max_bound, Hdeti));
  // const auto B = (ej - ei)*mKernel.kernelValue(min_bound, Hdeti);
  // return Vector(2.0*M_PI/(ei*ej)*FastMath::cube(Hdeti)*(A - B - Hdeti/ej*mInterp(min_bound, max_bound)));
}

}
