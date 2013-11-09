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
