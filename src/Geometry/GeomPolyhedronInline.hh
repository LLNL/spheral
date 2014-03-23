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
  return mVertexUnitNorms;
}

inline
const std::vector<std::vector<unsigned> >&
GeomPolyhedron::
vertexFacetConnectivity() const {
  return mVertexFacetConnectivity;
}

inline
const std::vector<std::vector<unsigned> >&
GeomPolyhedron::
facetFacetConnectivity() const {
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

}
