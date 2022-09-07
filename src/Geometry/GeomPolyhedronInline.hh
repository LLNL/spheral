#include "FacetedVolumeUtilities.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the vertices and facets.
//------------------------------------------------------------------------------
inline
const std::vector<GeomPolyhedron::Vector>&
GeomPolyhedron::
vertices() const {
  return mVertices;
}

inline
const std::vector<GeomPolyhedron::Facet>&
GeomPolyhedron::
facets() const {
  return mFacets;
}

inline
const std::vector<GeomPolyhedron::Vector>&
GeomPolyhedron::
vertexUnitNorms() const {
  if (mVertexUnitNorms.size() == 0) {
    GeometryUtilities::computeAncillaryGeometry(*this, 
                                                const_cast<std::vector<std::vector<unsigned>>&>(mVertexFacetConnectivity), 
                                                const_cast<std::vector<std::vector<unsigned>>&>(mFacetFacetConnectivity),
                                                const_cast<std::vector<Vector>&>(mVertexUnitNorms), 
                                                false);
  }
  return mVertexUnitNorms;
}

inline
const std::vector<std::vector<unsigned> >&
GeomPolyhedron::
vertexFacetConnectivity() const {
  if (mVertexFacetConnectivity.size() == 0) {
    GeometryUtilities::computeAncillaryGeometry(*this, 
                                                const_cast<std::vector<std::vector<unsigned>>&>(mVertexFacetConnectivity), 
                                                const_cast<std::vector<std::vector<unsigned>>&>(mFacetFacetConnectivity),
                                                const_cast<std::vector<Vector>&>(mVertexUnitNorms), 
                                                false);
  }
  return mVertexFacetConnectivity;
}

inline
const std::vector<std::vector<unsigned> >&
GeomPolyhedron::
facetFacetConnectivity() const {
  if (mFacetFacetConnectivity.size() == 0) {
    GeometryUtilities::computeAncillaryGeometry(*this, 
                                                const_cast<std::vector<std::vector<unsigned>>&>(mVertexFacetConnectivity), 
                                                const_cast<std::vector<std::vector<unsigned>>&>(mFacetFacetConnectivity),
                                                const_cast<std::vector<Vector>&>(mVertexUnitNorms), 
                                                true);
  }
  return mFacetFacetConnectivity;
}

//------------------------------------------------------------------------------
// Access the bounding box.
//------------------------------------------------------------------------------
inline
const GeomPolyhedron::Vector&
GeomPolyhedron::
xmin() const {
  return mXmin;
}

inline
const GeomPolyhedron::Vector&
GeomPolyhedron::
xmax() const {
  return mXmax;
}

//------------------------------------------------------------------------------
// Useful facet properties.
//------------------------------------------------------------------------------
inline
GeomPolyhedron::Vector
GeomPolyhedron::
facetCentroid(const unsigned facetID) const {
  REQUIRE(facetID < mFacets.size());
  return mFacets[facetID].position();
}

inline
double
GeomPolyhedron::
facetArea(const unsigned facetID) const {
  REQUIRE(facetID < mFacets.size());
  return mFacets[facetID].area();
}

inline
GeomPolyhedron::Vector
GeomPolyhedron::
facetAreaNormal(const unsigned facetID) const {
  REQUIRE(facetID < mFacets.size());
  return mFacets[facetID].area() * mFacets[facetID].normal();
}

}
