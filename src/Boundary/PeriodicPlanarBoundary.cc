//---------------------------------Spheral++----------------------------------//
// PeriodicPlanarBoundary -- The nested class member of Periodic Boundary that
// does the actual work.
//
// Created by JMO, Wed Apr 19 15:00:50 PDT 2000
//----------------------------------------------------------------------------//

// #include "PeriodicBoundary.hh"
#include "Utilities/DBC.hh"

using NodeSpace::NodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
PeriodicPlanarBoundary():
  PlanarBoundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given planes.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
PeriodicPlanarBoundary(const GeomPlane<Dimension>& plane1,
                       const GeomPlane<Dimension>& plane2):
  PlanarBoundary<Dimension>(plane1, plane2) {

  // Once we're done the boundary condition should be in a valid state.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
~PeriodicPlanarBoundary() {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to the ghost values in the given field.
// For periodic boundaries, this is always simple.  Just perform a copy of the
// control to the ghost values.
//------------------------------------------------------------------------------
// int fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, int>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 and *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Scalar fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 and *controlItr < nodeList.numNodes());
    CHECK2(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes(),
           "bounds error:  " << *ghostItr << " " << nodeList.firstGhostNode() << " " << nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Vector fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 and *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 and *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Symmetric Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 and *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

// Third Rank Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  vector<int>::const_iterator controlItr = this->controlBegin(nodeList);
  vector<int>::const_iterator ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr >= 0 and *controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes());
    field(*ghostItr) = field(*controlItr);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary condition on nodes in violation for the given field.
// For periodic boundaries, this is always simple.  Do nothing!
//------------------------------------------------------------------------------
// int fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, int>& field) const {
}

// Scalar fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Vector fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
}

// Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Symmetric Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

// Third Rank Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
}

//------------------------------------------------------------------------------
// Test if the periodic planar boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::valid() const {
  return (PlanarBoundary<Dimension>::valid());
}
