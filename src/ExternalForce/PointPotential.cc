//---------------------------------Spheral++----------------------------------//
// PointPotential -- Impose a potential from a point mass.
//
// Created by JMO, Sun Mar 30 22:08:55 PST 2003
//----------------------------------------------------------------------------//
#include "PointPotential.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PointPotential<Dimension>::
PointPotential(double G, double mass, double coreRadius, 
               const typename Dimension::Vector& origin):
  GenericBodyForce<Dimension>(),
  mG(G),
  mMass(mass),
  mCoreRadius2(coreRadius*coreRadius),
  mOrigin(origin),
  mDeltaPhiFraction(0.01),
  mPotentialEnergy(0.0) {
  ENSURE(mG > 0.0);
  ENSURE(mMass >= 0.0);
  ENSURE(mCoreRadius2 >= 0.0);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PointPotential<Dimension>::
~PointPotential() {
}

//------------------------------------------------------------------------------
// Calculate the acceleration due to the point mass on the given set of nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PointPotential<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the node positions from the state.
  const FieldList<Dimension, Scalar> mnode = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);

  // Get the acceleration and position change vectors we'll be modifying.
  FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Loop over the internal nodes.
  mPotentialEnergy = 0.0;
  const Scalar thpt = G()*mass();
  for (InternalNodeIterator<Dimension> nodeItr = dataBase.internalNodeBegin();
       nodeItr < dataBase.internalNodeEnd();
       ++nodeItr) {
    const Vector r = position(nodeItr) - mOrigin;
    const Vector runit = r.unitVector();
    const Scalar rsoft2 = r.magnitude2() + mCoreRadius2;
    const Scalar rsoft  = sqrt(rsoft2);
    const Scalar rsoft3 = rsoft2*rsoft;
    CHECK(rsoft3 > 0.0);
    CHECK(rsoft  > 0.0);
    DxDt(nodeItr) += velocity(nodeItr);
    DvDt(nodeItr) -= thpt*r.magnitude()/rsoft3*runit;
    mPotentialEnergy -= thpt*mnode(nodeItr)/rsoft;
  }
}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename PointPotential<Dimension>::TimeStepType
PointPotential<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const typename Dimension::Scalar currentTime) const {

  // This is a hard one to pick for this package.  We'll choose a timestep
  // such that no nodes potential should change more than a set fraction.
  // (Delta phi)/phi = -(Delta r)/(r^2 + rc^2)
  // dt = (Delta phi)/phi (r^2 + rc^2) / v
  Scalar mindt = FLT_MAX;
  Scalar minr, minv;
  FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  for (InternalNodeIterator<Dimension> nodeItr = dataBase.internalNodeBegin();
       nodeItr < dataBase.internalNodeEnd();
       ++nodeItr) {
    const Scalar rsoft = sqrt((position(nodeItr) - mOrigin).magnitude2() +
                              mCoreRadius2);
    const Scalar v = velocity(nodeItr).magnitude() + 1.0e-10;
    const Scalar dt = mDeltaPhiFraction*rsoft/v;
    if (dt < mindt) {
      mindt = dt;
      minr = rsoft;
      minv = v;
    }
  }

  std::stringstream reasonStream;
  reasonStream << "mindt = " << mindt << " | "
               << "rsoft = " << minr << " " 
               << "minv = " << minv << std::endl;

  return TimeStepType(mindt, reasonStream.str());
}

}
