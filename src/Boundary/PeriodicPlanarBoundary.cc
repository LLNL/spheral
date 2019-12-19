//---------------------------------Spheral++----------------------------------//
// PeriodicPlanarBoundary -- The nested class member of Periodic Boundary that
// does the actual work.
//
// Created by JMO, Wed Apr 19 15:00:50 PDT 2000
//----------------------------------------------------------------------------//

// #include "PeriodicBoundary.hh"
#include "Utilities/DBC.hh"

namespace {

//------------------------------------------------------------------------------
// Map a faceted volume
//------------------------------------------------------------------------------
template<typename BC, typename Dimension>
struct MapFacetedVolume {
  typename Dimension::FacetedVolume 
  static doit(const BC& bc,
       const typename Dimension::FacetedVolume& poly) {
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::FacetedVolume FacetedVolume;
    const auto& enterPlane = bc.enterPlane();
    const auto& exitPlane = bc.exitPlane();
    const auto& verts0 = poly.vertices();
    const auto& facets = poly.facetVertices();
    const auto n = verts0.size();
    vector<Vector> verts1;
    for (auto& vert: verts0) verts1.push_back(bc.mapPosition(vert, enterPlane, exitPlane));
    return FacetedVolume(verts1, facets);
  }
};

template<typename BC>
struct MapFacetedVolume<BC, Dim<1>> {
  Dim<1>::FacetedVolume 
  static doit(const BC& bc,
       const Dim<1>::FacetedVolume& poly) {
    const auto& enterPlane = bc.enterPlane();
    const auto& exitPlane = bc.exitPlane();
    return Dim<1>::FacetedVolume(bc.mapPosition(poly.center(),
                                                enterPlane,
                                                exitPlane),
                                 poly.extent());
  }
};

}

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

// Fourth Rank Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {

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

// Fifth Rank Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {

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

// FacetedVolume fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {

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
    field(*ghostItr) = MapFacetedVolume<PeriodicBoundary<Dimension>::PeriodicPlanarBoundary, Dimension>::doit(*this, field(*controlItr));
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary condition on nodes in violation for the given field.
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

// Fourth Rank Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
}

// Fifth Rank Tensor fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
}

// FacetedVolume fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  REQUIRE(valid());
  const auto& nodeList = field.nodeList();
  for (auto itr = this->violationBegin(nodeList);
       itr < this->violationEnd(nodeList); 
       ++itr) {
    CHECK(*itr >= 0 && *itr < nodeList.numInternalNodes());
    field(*itr) = MapFacetedVolume<PeriodicBoundary<Dimension>::PeriodicPlanarBoundary, Dimension>::doit(*this, field(*itr));
  }
}

//------------------------------------------------------------------------------
// Test if the periodic planar boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::valid() const {
  return (PlanarBoundary<Dimension>::valid());
}
