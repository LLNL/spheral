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
  double operator()(const Dim<2>::Vector& etas) const {
    const auto etaj = etas[0];
    const auto etai = etas[1];
    const auto low = std::abs(etaj - etai);
    const auto high = std::min(metaMax, etaj + etai);
    if (low >= high) return 0.0;
    return simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                        low,
                                                        high,
                                                        mW.numPoints());
  }
};

}

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
SphericalTableKernel::SphericalTableKernel(const TableKernel<Dim<3>>& kernel):
  mInterp(Dim<2>::Vector(0.0, 0.0),
          Dim<2>::Vector(5.0*kernel.kernelExtent(), 5.0*kernel.kernelExtent()),
          kernel.numPoints(),
          kernel.numPoints(),
          W3S1Func(kernel)),
  mGradInterp(),
  mGrad2Interp(),
  mKernel(kernel) {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
SphericalTableKernel::SphericalTableKernel(const SphericalTableKernel& rhs):
  mInterp(rhs.mInterp),
  mGradInterp(rhs.mGradInterp),
  mGrad2Interp(rhs.mGrad2Interp),
  mKernel(rhs.mKernel) {
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
    mGrad2Interp = rhs.mGrad2Interp;
    mKernel = rhs.mKernel;
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

}
