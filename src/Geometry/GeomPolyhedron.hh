//---------------------------------Spheral++----------------------------------//
// GeomPolyhedron -- Geometric polyhedron class.
//
// A 3-D structure representing a polyhedron as a collection of GeomFacets.
//
// Created by JMO, Fri Jan 29 14:36:01 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomPolyhedron__
#define __Spheral_GeomPolyhedron__

#include <vector>
#include "GeomVector.hh"
#include "GeomFacet3d.hh"

namespace Spheral {

class GeomPolyhedron {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef GeomVector<3> Vector;
  typedef GeomFacet3d Facet;

  //----------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //----------------------------------------------------------------------------
  GeomPolyhedron();

  // Note the following constructor constructs the convex hull of the given points,
  // meaning that the full set of points passed in may not appear in the vertices.
  GeomPolyhedron(const std::vector<Vector>& points);

  // This constructor takes a set of points that define facets of the surface
  // of the polyhedron.
  GeomPolyhedron(const std::vector<Vector>& points,
                 const std::vector<std::vector<unsigned> >& facetIndices);

  GeomPolyhedron(const GeomPolyhedron& rhs);
  GeomPolyhedron& operator=(const GeomPolyhedron& rhs);
  ~GeomPolyhedron();

  // Test if the given point is internal to the polyhedron.
  bool contains(const Vector& point,
                const bool countBoundary = true,
                const double tol = 1.0e-8) const;

  // If we have a convex polyhedron we can do a faster containment test.
  bool convexContains(const Vector& point,
                      const bool countBoundary = true,
                      const double tol = 1.0e-8) const;

  // Test if we intersect another polyhedron.
  bool intersect(const GeomPolyhedron& rhs) const;
  bool convexIntersect(const GeomPolyhedron& rhs) const;
  bool intersect(const std::pair<Vector, Vector>& rhs) const;  // Another way of representing a box.

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

  // Spit out a vector<vector<unsigned> > and vector<Vector> that encode the facets.
  std::vector<std::vector<unsigned> > facetVertices() const;
  std::vector<Vector> facetNormals() const;

  // Reconstruct the internal data given a set of verticies and the vertex
  // indicies that define the facets.
  void reconstruct(const std::vector<Vector>& vertices,
                   const std::vector<std::vector<unsigned> >& facetVertices,
                   const std::vector<Vector>& facetNormals);

  // Compute the volume.
  double volume() const;

  // Find the facet closest to the given point.
  unsigned closestFacet(const Vector& p) const;

  // Compute the minimum distance to a point.
  double distance(const Vector& p) const;

  // Find the point in the polyhedron closest to the given point.
  Vector closestPoint(const Vector& p) const;

  // Shift by a Vector delta.
  GeomPolyhedron& operator+=(const Vector& rhs);
  GeomPolyhedron& operator-=(const Vector& rhs);
  GeomPolyhedron operator+(const Vector& rhs) const;
  GeomPolyhedron operator-(const Vector& rhs) const;

  // Scale by a scalar.
  GeomPolyhedron& operator*=(const double rhs);
  GeomPolyhedron& operator/=(const double rhs);
  GeomPolyhedron operator*(const double rhs) const;
  GeomPolyhedron operator/(const double rhs) const;

  // Comparisons.
  bool operator==(const GeomPolyhedron& rhs) const;
  bool operator!=(const GeomPolyhedron& rhs) const;

  // Test if the polyhedron is convex.
  bool convex(const double tol = 1.0e-8) const;

  // Set the bounding box.
  void setBoundingBox();

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<Vector> mVertices;
  std::vector<Facet> mFacets;
  std::vector<Vector> mVertexUnitNorms;
  std::vector<std::vector<unsigned> > mVertexFacetConnectivity, mFacetFacetConnectivity;
  Vector mXmin, mXmax, mCentroid;
  double mRinterior2;
  bool mConvex;

  static FILE* mDevnull;
};

std::ostream& operator<<(std::ostream& os, const GeomPolyhedron& polygon);

}

#include "GeomPolyhedronInline.hh"

#endif

