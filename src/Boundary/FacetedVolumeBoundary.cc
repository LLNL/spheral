//---------------------------------Spheral++----------------------------------//
// FacetedVolumeBoundary -- Apply a FacetedVolume boundary condition to Spheral++
// Fields.
//----------------------------------------------------------------------------//
#include "findNodesTouchingThroughPlanes.hh"
#include "mapPositionThroughPlanes.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/innerProduct.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Mesh/Mesh.hh"

#include "FacetedVolumeBoundary.hh"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Construct the plane associated with a Facet
//------------------------------------------------------------------------------
GeomPlane<Dim<2>>
facetPlane(const GeomFacet2d& facet, const bool interiorBoundary) {
  const auto& p1 = facet.point1();
  const auto  nhat = (interiorBoundary ?
                       facet.normal() :
                      -facet.normal());
  return GeomPlane<Dim<2>>(p1, nhat);
}

GeomPlane<Dim<3>>
facetPlane(const GeomFacet3d& facet, const bool interiorBoundary) {
  const auto& p1 = facet.point(0);
  const auto  nhat = (interiorBoundary ?
                       facet.normal() :
                      -facet.normal());
  return GeomPlane<Dim<3>>(p1, nhat);
}

//------------------------------------------------------------------------------
// Find the node touching a given facet (Polygon)
//------------------------------------------------------------------------------
std::vector<int>
nodesTouchingFacet(const NodeList<Dim<2>>& nodes,
                   const Facet2d& facet,
                   const bool interiorBoundary) {
  const auto& p1 = facet.point1();
  const auto& p2 = facet.point2();
  const auto  phat = (p2 - p1).unitVector();
  const auto  smax = (p2 - p1).magnitude();
  const auto  plane = facetPlane(facet, interiorBoundary);
  const auto potentials = findNodesTouchingThroughPlanes(nodes, plane, plane);
  std::vector<result> result;
  const auto& pos = nodes.positions();
  for (const auto i: potentials) {
    const auto s = (pos(i) - p0).dot(phat);
    if (s >= 0.0 and s <= smax) result.push_back(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Find the node touching a given facet (Polyhedron)
//------------------------------------------------------------------------------
std::vector<int>
nodesTouchingFacet(const NodeList<Dim<3>>& nodes,
                   const Facet3d& facet,
                   const bool interiorBoundary) {
  std::vector<result> result;
  return result;
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FacetedVolumeBoundary<Dimension>::FacetedVolumeBoundary(const FacetedVolume& poly,
                                                        const bool interiorBoundary,
                                                        const bool useGhosts):
  Boundary<Dimension>(),
  mInteriorBoundary(interiorBoundary),
  mUseGhosts(useGhosts),
  mReflectOperators() {

  // Build the reflection operators for each facet of the poly.
  const auto& facets = poly.facets();
  for (const auto& facet: facets) {
    mReflectOperators.push_back(interiorBoundary ?
                                planarReflectingOperator(facet.normal()) :
                                planarReflectingOperator(-facet.normal()));
  }
  ENSURE(mReflectOperators.size() == facets.size());

}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FacetedVolumeBoundary<Dimension>::~FacetedVolumeBoundary() {
}

//------------------------------------------------------------------------------
// Create ghost nodes for the NodeList
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::setGhostNodes(NodeList<Dimension>& nodeList) {
  if (mUseGhosts) {
    // Remember which node list we are setting the ghost nodes for.
    const auto name = nodeList.name();
    this->addNodeList(nodeList);
    auto& boundaryNodes = this->accessBoundaryNodes(nodeList);
    mFacetControlNodes[nodeList.name()].clear();
    const auto& facets = mPoly.facets();
    for (const auto& facet: facets) {
      mFacetControlNodes[name].push_back(nodesTouchingFacet(nodeList, facet, mInteriorBoundary));
      boundaryNodes.controlNodes.insert(boundaryNodes.controlNodes.end(),
                                        mFacetControlNodes[name].back().begin(),
                                        mFacetControlNOdes[name].back().end());
    }
    CHECK(mFacetControlNodes[name].size() == facets.size());

    // Create the ghost nodes
    auto firstGhost = nodeList.numNodes();
    nodeList.numGhostNodes(nodeList.numGhostNodes() + boundaryNodes.controlNodes.size());
    boundaryNodes.ghostNodes.resize(boundaryNodes.controlNodes.size());
    mFacetGhostNodes[name].clear();
    for (const auto& controls: mFacetControlNodes[name]) {
      mFacetGhostNodes.push_back(std::make_pair(firstGhost, firstGhost + controls.size()));
      firstGhost += controls.size();
    }
    CHECK(mFacetGhostNodes[name].size() == facets.size());

    // Assign the positions and H's to the new ghost nodes.
    updateGhostNodes(nodeList);
  }
}

//------------------------------------------------------------------------------
// Internal method to set the minimal field values necessary to make the 
// ghost nodes valid.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::updateGhostNodes(NodeList<Dimension>& nodeList) {

  // Walk each of the facets
  const auto& facets = mPoly.facets();
  const auto  nfacets = facets.size();
  const auto  name = nodeList.name();
  const auto& controls = mFacetControlNodes[name];
  const auto& ghostRanges = mFacetGhostNodes[name];
  CHECK(controls.size() == nfacets);
  CHECK(ghostRanges.size() == nfacets);
  auto&       pos = nodeList.positions();
  for (auto ifacet = 0; ifacet < nfacets; ++ifacet) {
    const auto n = controls[ifacet].size();
    CHECK(ghostRanges[ifacet].second - ghostRanges[ifacet].first == n);

    for (auto k = 0; k < n; ++k) {
      const auto i = controls[ifacet][k];
      const auto j = ghostRanges[ifacet].first + k;
      

      pos(j) = mapPosition(pos(i),
                                       mExitPlane,
                                       mEnterPlane);
    // CHECK(positions(*ghostItr) <= mEnterPlane);
  }

  // Set the Hfield.
  Field<Dimension, SymTensor>& Hfield = nodeList.Hfield();
  this->applyGhostBoundary(Hfield);

  // // Update the neighbor information.
  // nodeList.neighbor().updateNodes(); // (ghostNodes);
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields, just perform a copy.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {
}

// Specialization for scalar fields, just perform a copy.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for Vector fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
}

// Specialization for Tensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

// Specialization for ThirdRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
}

// Specialization for FourthRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
}

// Specialization for FifthRankTensor fields.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
}

// Specialization for FacetedVolumes
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
}

//------------------------------------------------------------------------------
// Enforce the boundary condition on the set of nodes in violation of the 
// boundary.
//------------------------------------------------------------------------------
// Specialization for int fields.  A no-op.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
}

// Specialization for scalar fields.  A no-op.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for vector fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Specialization for tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

// Specialization for third rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
}

// Specialization for fourth rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
}

// Specialization for fifth rank tensor fields.  Apply the reflection operator.
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
}

// Specialization for FacetedVolumes
template<typename Dimension>
void
FacetedVolumeBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
}

}
