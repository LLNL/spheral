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
  return mVertexUnitNorms;
}

inline
const std::vector<std::vector<unsigned> >&
GeomPolygon::
vertexFacetConnectivity() const {
  return mVertexFacetConnectivity;
}

inline
const std::vector<std::vector<unsigned> >&
GeomPolygon::
facetFacetConnectivity() const {
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

}
