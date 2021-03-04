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
    double operator()(const double eta) const {
      return eta * mW(eta, 1.0);
    }
  };

  // Now the lookup based on (etaj, etai) -- the (rprime/h, r/h) from the paper
  double operator()(const Dim<2>::Vector& bounds) const {
    const auto low = bounds[0];
    const auto high = bounds[1];
    if (low >= high) return 0.0;
    return simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                        low,
                                                        high,
                                                        std::max(size_t(1000), mW.numPoints()));
  }
};

//------------------------------------------------------------------------------
// Functor for doing our volume integration of the kernel gradient
//------------------------------------------------------------------------------
struct gradW3S1Func {
  const TableKernel<Dim<3>>& mW;
  double metaMax;
  gradW3S1Func(const TableKernel<Dim<3>>& W): mW(W), metaMax(mW.kernelExtent()) {}

  // Define a local nested functor we'll use to do the volume integral for
  // an (etaj, etai) lookup.
  struct VolFunc {
    const TableKernel<Dim<3>>& mW;
    VolFunc(const TableKernel<Dim<3>>& W): mW(W) {}
    double operator()(const double eta) const {
      return eta * mW.gradValue(eta, 1.0) + mW(eta, 1.0);
    }
  };

  // Now the lookup based on (etaj, etai) -- the (rprime/h, r/h) from the paper
  double operator()(const Dim<2>::Vector& bounds) const {
    const auto low = bounds[0];
    const auto high = bounds[1];
    if (low >= high) return 0.0;
    return simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                        low,
                                                        high,
                                                        std::max(size_t(1000), mW.numPoints()));
  }
};

}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
SphericalTableKernel::SphericalTableKernel(const TableKernel<Dim<3>>& kernel):
  mInterp(Dim<2>::Vector(0.0, 0.0),
          Dim<2>::Vector(kernel.kernelExtent(), kernel.kernelExtent()),
          std::max(size_t(200), kernel.numPoints()),
          std::max(size_t(200), kernel.numPoints()),
          W3S1Func(kernel)),
  // mGradInterp(Dim<2>::Vector(0.0, 0.0),
  //             Dim<2>::Vector(kernel.kernelExtent(), kernel.kernelExtent()),
  //             std::max(size_t(200), kernel.numPoints()),
  //             std::max(size_t(200), kernel.numPoints()),
  //             gradW3S1Func(kernel)),
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
//------------------------------------------------------------------------------
double
SphericalTableKernel::grad(const Dim<1>::Vector& etaj,
                           const Dim<1>::Vector& etai,
                           const Dim<1>::Scalar  Hdeti) const {
  REQUIRE(Hdeti >= 0.0);
  const auto ei = std::max(1e-10, etai[0]);
  const auto ej = std::max(1e-10, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) return 0.0;
  const auto max_bound = std::min(metamax, ei + ej);
  const auto etahat = sgn(ei - ej);
  const auto A = (ei + ej >= metamax ?
                  0.0 :
                  max_bound*mKernel.kernelValue(max_bound, Hdeti));
  const auto B = min_bound*mKernel.kernelValue(min_bound, Hdeti)*etahat;
  return 2.0*M_PI/(ei*ej)*Hdeti*(A + B - Hdeti/ei*mInterp(Dim<2>::Vector(min_bound, max_bound)));

  // const auto A = (ei + ej < metamax ?
  //                 Hdeti*(max_bound*mKernel.gradValue(max_bound, Hdeti) + mKernel.kernelValue(max_bound, Hdeti)) :
  //                 0.0);
  // const auto B = Hdeti*(ei - ej)/std::max(1.0e-10, min_bound)*(min_bound*mKernel.gradValue(min_bound, Hdeti) + mKernel.kernelValue(min_bound, Hdeti));
  // // const auto B = Hdeti*(ei - ej)/std::max(1.0e-10, min_bound)*(min_bound*mKernel.gradValue(min_bound, Hdeti) + mKernel.kernelValue(min_bound, Hdeti));
  // // const auto B = (min_bound*Hdeti*mKernel.gradValue(min_bound, Hdeti) + Hdeti*mKernel.kernelValue(min_bound, Hdeti) * (ei >= ej ? 1.0 : -1.0);
  // return 2.0*M_PI/(ei*ej)*Hdeti*((A - B) - Hdeti*mInterp(Dim<2>::Vector(min_bound, max_bound)/ei));

  // // const auto grad_min_bound = ej > ei ?           1.0 : -1.0;
  // // const auto grad_max_bound = ei + ej < metamax ? 1.0 :  0.0;
  // // return 2.0*M_PI/(ei*ej)*Hdeti*(-mInterp(Dim<2>::Vector(min_bound, max_bound))*Hdeti/ej +
  // //                                Hdeti*((max_bound*mKernel.gradValue(max_bound, Hdeti) + mKernel.kernelValue(max_bound, Hdeti))*grad_max_bound -
  // //                                       (min_bound*mKernel.gradValue(min_bound, Hdeti) + mKernel.kernelValue(min_bound, Hdeti))*grad_min_bound));
  // // return 2.0*M_PI*Hdeti*Hdeti*(mGradInterp(Dim<2>::Vector(min_bound, max_bound))/(ei*ej)*sgn(ei - ej) -
  // //                              mInterp(Dim<2>::Vector(min_bound, max_bound))*Hdeti/(ei*ej*ej));
}

}
