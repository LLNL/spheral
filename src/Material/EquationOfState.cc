//---------------------------------Spheral++----------------------------------//
// EquationOfState -- Abstract base class for the equation of state classes.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//

#include "EquationOfState.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/bisectRoot.hh"

namespace Spheral {
namespace Material {

namespace {

template<typename Dimension>
struct Pfunctor {
  const EquationOfState<Dimension>& mEOS;
  NodeSpace::NodeList<Dimension> mNodes;
  mutable FieldSpace::Field<Dimension, double> mRho, mEps, mP;

  Pfunctor(const EquationOfState<Dimension>& eos, const double rho):
    mEOS(eos),
    mNodes("__dummynodes__", 1, 0),
    mRho("rho", mNodes, rho),
    mEps("eps", mNodes),
    mP("P", mNodes) {}

  double operator()(const double eps) const {
    mEps[0] = eps;
    mEOS.setPressure(mP, mRho, mEps);
    return mP[0];
  }
};

}

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
EquationOfState<Dimension>::EquationOfState(const PhysicalConstants& constants,
                                            const double minimumPressure,
                                            const double maximumPressure,
                                            const MaterialPressureMinType minPressureType):
  mConstants(constants),
  mMinimumPressure(minimumPressure),
  mMaximumPressure(maximumPressure),
  mMinPressureType(minPressureType) {
  REQUIRE(mMinimumPressure <= mMaximumPressure);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
EquationOfState<Dimension>::~EquationOfState() {
}

//------------------------------------------------------------------------------
// Look up an energy that gives the requested pressure at the specified density.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
EquationOfState<Dimension>::
specificThermalEnergyForPressure(const typename Dimension::Scalar Ptarget,
                                 const typename Dimension::Scalar rho,
                                 const typename Dimension::Scalar epsMin,
                                 const typename Dimension::Scalar epsMax,
                                 const typename Dimension::Scalar epsTol,
                                 const typename Dimension::Scalar Ptol,
                                 const unsigned maxIterations) const {
  const Pfunctor<Dimension> pfunc(*this, rho);
  return bisectRoot(pfunc, epsMin, epsMax, epsTol, Ptol, maxIterations);
}

}
}
