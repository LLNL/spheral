//---------------------------------Spheral++----------------------------------//
// GeomFacet2d -- A facet of a polygon (just two points).
//
// Note a Facet does not maintain it's own copies of it's end points -- the
// assumption is that this is a Facet of a GeomPolygon and that polygon owns
// the set of vertex positions.
//
// Created by JMO, Thu Jan 28 10:58:32 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomFacet2d__
#define __Spheral_GeomFacet2d__

#include <vector>
#include <iostream>

#include "Geometry/GeomVector_fwd.hh"

namespace Spheral {

// Forward declare the polygon.
class GeomPolygon;

class GeomFacet2d {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef GeomVector<2> Vector;

  //----------------------------------------------------------------------------
  // Constructors, assignment, destructor.
  //----------------------------------------------------------------------------
  GeomFacet2d();
  GeomFacet2d(const std::vector<Vector>& vertices,
              const unsigned point1,
              const unsigned point2);
  GeomFacet2d(const GeomFacet2d& rhs);
  GeomFacet2d& operator=(const GeomFacet2d& rhs);
  ~GeomFacet2d();

  // Is the given point above, below, or colinear with the facet?
  int compare(const Vector& point,
              const double tol = 1.0e-8) const;

  // Compare a set of points:
  //  1 => all points above.
  //  0 => points both above and below (or equal).
  // -1 => all points below.
  int compare(const std::vector<Vector>& points,
              const double tol = 1.0e-8) const;

  // Access the points.
  const Vector& point1() const;
  const Vector& point2() const;

  unsigned ipoint1() const;
  unsigned ipoint2() const;

  // We provide the points as a standard vector for compatibility with GeomFacet3d.
  const std::vector<unsigned>& ipoints() const;

  const Vector& normal() const;
  
  Vector position() const;

  double area() const;

  // Compute the minimum distance from the facet to a point.
  double distance(const Vector& p) const;

  // Compute the closest point on the facet to the given point.
  Vector closestPoint(const Vector& p) const;

  // Comparisons.
  bool operator==(const GeomFacet2d& rhs) const;
  bool operator!=(const GeomFacet2d& rhs) const;

  // For consistency with GeomFacet3d.
  void decompose(std::vector<std::array<Vector, 2>>& subfacets) const;
  
private:
  //--------------------------- Private Interface ---------------------------//
  const std::vector<Vector>* mVerticesPtr;
  std::vector<unsigned> mPoints;
  Vector mNormal;

  // Make polygon a friend.
  friend GeomPolygon;
};

// Provide an ostream operator for GeomFacet2d.
std::ostream& operator <<(std::ostream& os, const GeomFacet2d& facet);

}

#include "GeomFacet2dInline.hh"

#else 

namespace Spheral {
  class GeomFacet2d;
}

#endif

