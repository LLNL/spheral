//---------------------------------Spheral++----------------------------------//
// ThreeOneDimensionBoundary
// A cheezy way to run 1-D simulations using 3-D objects.  This simply ensures
// that Vectors have 0 x & y components, while tensor are diagonal with the
// yy and zz components set equal to the xx value.
//
// Created by JMO, Wed Jul 28 13:38:28 2004
//----------------------------------------------------------------------------//

#include "ThreeOneDimensionBoundary.hh"
#include "NodeList/FluidNodeList.hh"

#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ThreeOneDimensionBoundary<Dimension>::ThreeOneDimensionBoundary():
  Boundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ThreeOneDimensionBoundary<Dimension>::~ThreeOneDimensionBoundary() {
}

//------------------------------------------------------------------------------
// Determine the set of ghost nodes for the boundary condition.
// This boundary never has any ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::setGhostNodes(NodeList<Dimension>& nodeList) {
  // Add this NodeList, creating space for control & ghost nodes.
  addNodeList(nodeList);
}

//------------------------------------------------------------------------------
// Update the ghost nodes.  Also a no-op.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::updateGhostNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// Find the set of nodes in the given NodeList that violate this boundary
// condition.  Technically we flag all internal nodes as "violation" nodes,
// since we want to make sure they all meet our fake "1-D" approximation.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::setViolationNodes(NodeList<Dimension>& nodeList) {

  // Get the BoundaryNodes.violationNodes for this NodeList.
  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;
  addNodeList(nodeList);
  BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
  vector<int>& vNodes = boundaryNodes.violationNodes;
  vNodes.resize(0);

  // Set all internal nodes as violation nodes.
  for (int i = 0; i != nodeList.numInternalNodes(); ++i) vNodes.push_back(i);
  CHECK(vNodes.size() == nodeList.numInternalNodes());

  // Update the positions and H for the nodes in violation.
  updateViolationNodes(nodeList);
}

//------------------------------------------------------------------------------
// Update the nodes in violation of the boundary condition, fixing their
// positions and H tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::updateViolationNodes(NodeList<Dimension>& nodeList) {

  // Zero the y,z position values.
  Field<Dimension, Vector>& positions = nodeList.positions();
  this->enforceBoundary(positions);

  // Force the H tensors to be diagonal, with yy = zz = 1.
  Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
  for (int i = 0; i != Hfield.numElements(); ++i) {
    const double xx = Hfield(i).xx();
    Hfield(i).Zero();
    Hfield(i).xx(xx);
    for (int j = 1; j < Dimension::nDim; ++j) Hfield(i)(j,j) = 1.0;
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
// These are all no-ops for this boundary condition.
//------------------------------------------------------------------------------
// Specialization for scalar fields, just perform a copy.
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for Vector fields.
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
}

// Specialization for Tensor fields.
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

//------------------------------------------------------------------------------
// Enforce the boundary condition on the set of nodes in violation of the 
// boundary.
//------------------------------------------------------------------------------
// Specialization for scalar fields.  A no-op.
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for vector fields.  Zero out the x & y components.
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (int i = 0; i != field.numElements(); ++i) {
    const double x = field(i).x();
    field(i).Zero();
    field(i).x(x);
  }
}

// Specialization for tensor fields.  Zero the off diagonal terms, and set the yy
// & zz elements equal to xx.
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (int i = 0; i != field.numElements(); ++i) {
    const double xx = field(i).xx();
    field(i).Zero();
    for (int j = 0; j != Dimension::nDim; ++j) field(i)(j,j) = xx;
  }
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
ThreeOneDimensionBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {

  const NodeList<Dimension>& nodeList = field.nodeList();
  for (int i = 0; i != field.numElements(); ++i) {
    const double xx = field(i).xx();
    field(i).Zero();
    for (int j = 0; j != Dimension::nDim; ++j) field(i)(j,j) = xx;
  }
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

namespace Spheral {
template class ThreeOneDimensionBoundary< Dim<3> >;
}
