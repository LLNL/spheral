//---------------------------------Spheral++----------------------------------//
// GeomPolygon -- Geometric polygon class.
//
// A 2-D structure representing a polygon as a collection of GeomFacets.
//
// Created by JMO, Thu Jan 28 11:03:27 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomPolygon__
#define __Spheral_GeomPolygon__

#include <vector>

#include "GeomVector.hh"
#include "GeomFacet2d.hh"

namespace Spheral {

class GeomPolygon {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef GeomVector<2> Vector;
  typedef GeomFacet2d Facet;

  //----------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //----------------------------------------------------------------------------
  GeomPolygon();

  // Note the following constructor constructs the convex hull of the given points,
  // meaning that the full set of points passed in may not appear in the vertices.
  GeomPolygon(const std::vector<Vector>& points);

  // This constructor takes a set of points that define facets of the surface
  // of the polygon.
  GeomPolygon(const std::vector<Vector>& points,
              const std::vector<std::vector<unsigned> >& facetIndices);

  GeomPolygon(const GeomPolygon& rhs);
  GeomPolygon& operator=(const GeomPolygon& rhs);
  ~GeomPolygon();

  // Test if the given point is internal to the polygon.
  // This method works for any polygon.
  bool contains(const Vector& point,
                const bool countBoundary = true,
                const double tol = 1.0e-8) const;

  // If we have a convex polygon we can do a faster containment test.
  bool convexContains(const Vector& point,
                      const bool countBoundary = true,
                      const double tol = 1.0e-8) const;

  // Test if we intersect another polygon.
  bool intersect(const GeomPolygon& rhs) const;
  bool convexIntersect(const GeomPolygon& rhs) const;
  bool intersect(const std::pair<Vector, Vector>& rhs) const;  // Another way of representing a box.

  // Return the intersections of this polygon with a line segment 
  // denoted by it's end points.
  std::vector<Vector> intersect(const Vector& x0, const Vector& x1) const;

  // Return the centroid of the vertices.
  Vector centroid() const;

  // Access the internal data.
  const std::vector<Vector>& vertices() const;
  const std::vector<Facet>& facets() const;
  const std::vector<Vector>& vertexUnitNorms() const;
  const std::vector<std::vector<unsigned> >& vertexFacetConnectivity() const;
  const std::vector<std::vector<unsigned> >& facetFacetConnectivity() const;
  const Vector& xmin() const;
  const Vector& xmax() const;

  // Get the edges as integer (node) pairs.
  std::vector<std::pair<unsigned, unsigned> > edges() const;

  // Spit out a vector<vector<unsigned> > that encodes the facets.
  std::vector<std::vector<unsigned> > facetVertices() const;

  // Reconstruct the internal data given a set of verticies and the vertex
  // indicies that define the facets.
  void reconstruct(const std::vector<Vector>& vertices,
                   const std::vector<std::vector<unsigned> >& facetVertices);

  // Compute the volume.
  double volume() const;

  // Find the facet closest to the given point.
  unsigned closestFacet(const Vector& p) const;

  // Compute the minimum distance to a point.
  double distance(const Vector& p) const;

  // Find the point in the polyhedron closest to the given point.
  Vector closestPoint(const Vector& p) const;

  // Shift by a Vector delta.
  GeomPolygon& operator+=(const Vector& rhs);
  GeomPolygon& operator-=(const Vector& rhs);
  GeomPolygon operator+(const Vector& rhs) const;
  GeomPolygon operator-(const Vector& rhs) const;

  // Scale by a scalar.
  GeomPolygon& operator*=(const double rhs);
  GeomPolygon& operator/=(const double rhs);
  GeomPolygon operator*(const double rhs) const;
  GeomPolygon operator/(const double rhs) const;
  // Comparisons.
  bool operator==(const GeomPolygon& rhs) const;
  bool operator!=(const GeomPolygon& rhs) const;

  // Test if the polygon is convex.
  bool convex(const double tol = 1.0e-8) const;

  // Set the bounding box.
  void setBoundingBox();

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<Vector> mVertices;
  std::vector<Facet> mFacets;
  std::vector<Vector> mVertexUnitNorms;
  std::vector<std::vector<unsigned> > mVertexFacetConnectivity, mFacetFacetConnectivity;
  Vector mXmin, mXmax;
  bool mConvex;

  static FILE* mDevnull;
};

std::ostream& operator<<(std::ostream& os, const GeomPolygon& polygon);

}

#include "GeomPolygonInline.hh"

#endif

