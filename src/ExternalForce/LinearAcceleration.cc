//---------------------------------Spheral++----------------------------------//
// LinearAcceleration -- Impose a linear acceleration as a function of position
// on a given set of nodes.
//
// Created by JMO, Wed Sep 29 15:45:23 2004
//----------------------------------------------------------------------------//
#include "LinearAcceleration.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "DBC.hh"

namespace Spheral {
namespace PhysicsSpace {

using namespace std;
using NodeSpace::NodeList;
using FieldSpace::Field;
using DataBaseSpace::DataBase;
using FieldSpace::FieldList;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
LinearAcceleration<Dimension>::
LinearAcceleration(const Scalar a0,
                   const Scalar aslope):
  GenericBodyForce<Dimension>(),
  ma0(a0),
  maslope(aslope) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
LinearAcceleration<Dimension>::
~LinearAcceleration() {
}

//------------------------------------------------------------------------------
// Apply the acceleration to the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearAcceleration<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the state.
  const FieldList<Dimension, Vector> r = state.fields(HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Increment the acceleration.
  for (InternalNodeIterator<Dimension> itr = DvDt.internalNodeBegin();
       itr != DvDt.internalNodeEnd();
       ++itr) {
    const Scalar x = r(itr).x();
    Vector da;
    da.x(ma0 + maslope*x);
    DvDt(itr) += da;
  }
}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename LinearAcceleration<Dimension>::TimeStepType
LinearAcceleration<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const typename Dimension::Scalar currentTime) const {

  return TimeStepType(FLT_MAX, "No vote.");
}

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace PhysicsSpace {
template class LinearAcceleration< Dim<1> >;
template class LinearAcceleration< Dim<2> >;
template class LinearAcceleration< Dim<3> >;
}
}
