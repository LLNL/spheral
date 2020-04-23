//---------------------------------Spheral++----------------------------------//
// GeomFacet3d -- A facet of a polyhedron (triangular)
//
// Note a Facet does not maintain it's own copies of its vertices -- the
// assumption is that this is a Facet of a GeomPolyhedron and that polyhedron
// owns the set of vertex positions.
//
// Created by JMO, Fri Jan 29 14:29:28 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomFacet3d__
#define __Spheral_GeomFacet3d__

#include <vector>
#include <iostream>

#include "Geometry/GeomVector_fwd.hh"

namespace Spheral {

// Forward declare the polyhedron.
class GeomPolyhedron;  

class GeomFacet3d {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef GeomVector<3> Vector;

  //----------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //----------------------------------------------------------------------------
  GeomFacet3d();
  GeomFacet3d(const std::vector<Vector>& vertices,
              const std::vector<unsigned>& ipoints,
              const Vector& normal);
  GeomFacet3d(const GeomFacet3d& rhs);
  GeomFacet3d& operator=(const GeomFacet3d& rhs);
  ~GeomFacet3d();

  // Is the given point above, below, or coplanar with the facet?
  int compare(const Vector& point,
              const double tol = 1.0e-8) const;

  // Compare a set of points:
  //  1 => all points above.
  //  0 => points both above and below (or equal).
  // -1 => all points below.
  int compare(const std::vector<Vector>& points,
              const double tol = 1.0e-8) const;

  // Access the points and normal.
  const Vector& point(const unsigned index) const;
  const std::vector<unsigned>& ipoints() const;
  const Vector& normal() const;
  
  // Compute the position.
  Vector position() const;

  // Compute the area of the facet.
  double area() const;

  // Compute the minimum distance from the facet to a point.
  double distance(const Vector& p) const;

  // Compute the closest point on the facet to the given point.
  Vector closestPoint(const Vector& p) const;

  // Comparisons.
  bool operator==(const GeomFacet3d& rhs) const;
  bool operator!=(const GeomFacet3d& rhs) const;

  // Decompose the facet into triangles.
  void decompose(std::vector<std::array<Vector, 3>>& subfacets) const;
  
  // Split into triangular sub-facets.
  std::vector<GeomFacet3d> triangles() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const std::vector<Vector>* mVerticesPtr;
  std::vector<unsigned> mPoints;
  Vector mNormal;

  // Make Polyhedron a friend.
  friend GeomPolyhedron;
};

// Provide an ostream operator for GeomFacet3d.
std::ostream& operator <<(std::ostream& os, const GeomFacet3d& facet);

}

#include "GeomFacet3dInline.hh"

#else

namespace Spheral {
  class GeomFacet3d;
}

#endif

