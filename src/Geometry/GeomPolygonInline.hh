#include "FacetedVolumeUtilities.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the vertices and facets.
//------------------------------------------------------------------------------
inline
const std::vector<GeomPolygon::Vector>&
GeomPolygon::
vertices() const {
  return mVertices;
}

inline
const std::vector<GeomPolygon::Facet>&
GeomPolygon::
facets() const {
  return mFacets;
}

inline
const std::vector<GeomPolygon::Vector>&
GeomPolygon::
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
GeomPolygon::
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
GeomPolygon::
facetFacetConnectivity() const {
  if (mFacetFacetConnectivity.size() == 0) {  // lazy evaluation
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
const GeomPolygon::Vector&
GeomPolygon::
xmin() const {
  return mXmin;
}

inline
const GeomPolygon::Vector&
GeomPolygon::
xmax() const {
  return mXmax;
}

//------------------------------------------------------------------------------
// Useful facet properties.
//------------------------------------------------------------------------------
inline
double
GeomPolygon::
facetArea(const unsigned facetID) const {
  REQUIRE(facetID < mFacets.size());
  return mFacets[facetID].area();
}

inline
GeomPolygon::Vector
GeomPolygon::
facetAreaNormal(const unsigned facetID) const {
  REQUIRE(facetID < mFacets.size());
  return mFacets[facetID].area() * mFacets[facetID].normal();
}

}
