//---------------------------------Spheral++----------------------------------//
// ConstantAcceleration -- Impose a constant acceleration on a given set of
// nodes.
//
// Created by JMO, Mon Sep 27 23:01:15 PDT 2004
//----------------------------------------------------------------------------//
#include "ConstantAcceleration.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace PhysicsSpace {

using namespace std;
using NodeSpace::NodeList;
using FieldSpace::Field;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantAcceleration<Dimension>::
ConstantAcceleration(const Vector a0,
                     const NodeList<Dimension>& nodeList,
                     const vector<int>& indices):
  GenericBodyForce<Dimension>(),
  ma0(a0),
  mNodeListPtr(&nodeList),
  mIndices(indices) {
  for (vector<int>::const_iterator itr = mIndices.begin();
       itr != mIndices.end();
       ++itr) ENSURE(*itr >= 0 && *itr < mNodeListPtr->numNodes());
}

//------------------------------------------------------------------------------
// Constructor for all nodes in NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantAcceleration<Dimension>::
ConstantAcceleration(const Vector a0,
                     const NodeList<Dimension>& nodeList):
  GenericBodyForce<Dimension>(),
  ma0(a0),
  mNodeListPtr(&nodeList),
  mIndices() {
  mIndices.reserve(nodeList.numInternalNodes());
  for (unsigned i = 0; i != nodeList.numInternalNodes(); ++i) mIndices.push_back(i);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantAcceleration<Dimension>::
~ConstantAcceleration() {
}

//------------------------------------------------------------------------------
// Apply the acceleration to the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantAcceleration<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the acceleration state we'll be modifying.
  const typename State<Dimension>::KeyType key = State<Dimension>::buildFieldKey(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity,
                                                                                 mNodeListPtr->name());
  Field<Dimension, Vector>& DvDt = derivs.field(key, Vector::zero);

  // Increment the acceleration.
  for (vector<int>::const_iterator itr = mIndices.begin();
       itr != mIndices.end();
       ++itr) {
    CHECK(*itr >= 0 && *itr < mNodeListPtr->numNodes());
    DvDt(*itr) += ma0;
  }
}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename ConstantAcceleration<Dimension>::TimeStepType
ConstantAcceleration<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const typename Dimension::Scalar currentTime) const {

  return TimeStepType(FLT_MAX, "No vote.");
}

}
}

