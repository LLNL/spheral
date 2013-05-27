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
