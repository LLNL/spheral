//---------------------------------Spheral++----------------------------------//
// SphericalKernel
//
// Take a 3D Kernel and build a specialized 1D tabulated version appropriate
// for use with the spherical SPH algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Created by JMO, Wed Dec  2 16:41:20 PST 2020
//----------------------------------------------------------------------------//

#include "SphericalKernel.hh"
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

}            // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
SphericalKernel::SphericalKernel(const TableKernel<Dim<3>>& kernel):
  mInterp(0.0, kernel.kernelExtent(),
          0.0, kernel.kernelExtent(),
          std::max(size_t(200), kernel.numPoints()),
          std::max(size_t(200), kernel.numPoints()),
          W3S1Func(kernel)),
  mBaseKernel3d(kernel),
  mBaseKernel1d(kernel, kernel.numPoints()),
  metamax(kernel.kernelExtent()) {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
SphericalKernel::SphericalKernel(const SphericalKernel& rhs):
  mInterp(rhs.mInterp),
  mBaseKernel3d(rhs.mBaseKernel3d),
  mBaseKernel1d(rhs.mBaseKernel1d),
  metamax(rhs.metamax) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalKernel::~SphericalKernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
SphericalKernel&
SphericalKernel::operator=(const SphericalKernel& rhs) {
  if (this != &rhs) {
    mInterp = rhs.mInterp;
    mBaseKernel3d = rhs.mBaseKernel3d;
    mBaseKernel1d = rhs.mBaseKernel1d;
    metamax = rhs.metamax;
  }
  return *this;
}

}

// We need to instantiate the special TableKernel constructor we use
#include "TableKernel.cc"
namespace Spheral {
template TableKernel<Dim<1>>::TableKernel(const TableKernel<Dim<3>>&, const unsigned, const double);
}

