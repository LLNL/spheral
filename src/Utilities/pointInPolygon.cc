//------------------------------------------------------------------------------
// pointInPolygon
//
// Test if a given point is in a polygon or not.  This is based on code I found
// at
//
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//
// which I have modified here somewhat to use Spheral structures and such.
//------------------------------------------------------------------------------
#include "pointInPolygon.hh"
#include "testBoxIntersection.hh"
#include "pointOnPolygon.hh"
#include "pointDistances.hh"
#include "lineSegmentIntersections.hh"

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// Work with a closed polygon defined by it's vertices.
// We assume the caller has ordered the vertices (clockwise or 
// counter-clockwise), and should have already performed the box exclusion
// test for efficiency.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<2>::Vector& p,
                    const std::vector<Dim<2>::Vector>& vertices,
                    const bool countBoundary,
                    const double tol) {
  typedef Dim<2>::Vector Vector;

  // Check if the point is on the boundary (within tolerance).
  // The succeeding code sometimes (but uniquely) includes boundary points, so 
  // we need to check for boundary first.
  if (pointOnPolygon(p, vertices, tol)) return countBoundary;

  // Now we do the test of casting a semi-infinite ray in the x direction from 
  // the point and counting intersections with the polygon.
  // Note this bit of code will reject points lying exactly on the boundary
  // of a polygon.
  unsigned i, j;
  bool result = false;
  const unsigned nvertices = vertices.size();
  const double px = p.x();
  const double py = p.y();
  for (i = 0, j = nvertices - 1; i < nvertices; j = i++) {
    if ( ((vertices[i].y() > py) != (vertices[j].y() > py)) &&
	 (px < (vertices[j].x() - vertices[i].x()) * (py - vertices[i].y()) / (vertices[j].y() - vertices[i].y()) + vertices[i].x()) )
      result = not result;
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Test a polygon.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<2>::Vector& p,
                    const Dim<2>::FacetedVolume& polygon,
                    const bool countBoundary,
                    const double tol) {
  if (not testPointInBox(p, polygon.xmin(), polygon.xmax(), tol)) return false;
  return pointInPolygon(p, polygon.vertices(), countBoundary, tol);
}

//------------------------------------------------------------------------------
// Test a polygon in 3-D.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<3>::Vector& p,
                    const vector<Dim<3>::Vector>& vertices,
                    const Dim<3>::Vector& normal) {
  typedef Dim<3>::Vector Vector;

  // Prerequisites.
  const unsigned npts = vertices.size();
  unsigned i, j, k;
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(vertices.size() > 2);
    const double normmag = normal.magnitude();
    Vector normi;
    for (i = 1; i != npts; ++i) {
      j = (i + 1) % npts;
      k = (i + 2) % npts;
      normi = (vertices[j] - vertices[i]).cross(vertices[k] - vertices[i]);
      REQUIRE(fuzzyEqual(abs(normi.dot(normal)), normmag*normi.magnitude(), 1.0e-10));
    }
    REQUIRE(fuzzyEqual(pointPlaneDistance(p, vertices[0], normal.unitVector()), 0.0, 1.0e-10));
  }
  END_CONTRACT_SCOPE;

  bool result = false;

  const double px = p.x();
  const double py = p.y();
  const double pz = p.z();

  double fxmin =  1e50, fymin =  1e50, fzmin =  1e50;
  double fxmax = -1e50, fymax = -1e50, fzmax = -1e50;
  for (i = 0; i != npts; ++i) {
    fxmin = min(fxmin, vertices[i].x());
    fymin = min(fymin, vertices[i].y());
    fzmin = min(fzmin, vertices[i].z());
    fxmax = max(fxmax, vertices[i].x());
    fymax = max(fymax, vertices[i].y());
    fzmax = max(fzmax, vertices[i].z());
  }
  if (px >= fxmin and px <= fxmax and
      py >= fymin and py <= fymax and
      pz >= fzmin and pz <= fzmax) {

    // Check if the point is on the boundary.
    for (i = 0; i != npts; ++i) {
      j = (i + 1) % npts;
      if ((closestPointOnSegment(p, vertices[i], vertices[j]) - p).magnitude2() < 1.0e-10) return true;
    }

    if (abs(normal.x()) >= abs(normal.y()) and abs(normal.x()) >= abs(normal.z())) {

      // x plane -- use (y,z) coordinates.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[i].z() > pz) != (vertices[j].z() > pz)) &&
             (py < (vertices[j].y() - vertices[i].y()) * (pz - vertices[i].z()) / (vertices[j].z() - vertices[i].z()) + vertices[i].y()) )
          result = not result;
      }

    } else if (abs(normal.y()) >= abs(normal.x()) and abs(normal.y()) >= abs(normal.z())) {

      // y plane -- use (x,z) coordinates.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[i].z() > pz) != (vertices[j].z() > pz)) &&
             (px < (vertices[j].x() - vertices[i].x()) * (pz - vertices[i].z()) / (vertices[j].z() - vertices[i].z()) + vertices[i].x()) )
          result = not result;
      }

    } else {
      CHECK(abs(normal.z()) >= abs(normal.x()) and abs(normal.z()) >= abs(normal.y()));

      // z plane -- use (x,y) coordinate.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[i].y() > py) != (vertices[j].y() > py)) &&
             (px < (vertices[j].x() - vertices[i].x()) * (py - vertices[i].y()) / (vertices[j].y() - vertices[i].y()) + vertices[i].x()) )
          result = not result;
      }

    }
  }

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Test a polygon in 3-D.
// This version allows you to pass in an indirect addressing definition of the
// polygon vertices.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<3>::Vector& p,
                    const vector<Dim<3>::Vector>& vertices,
                    const vector<unsigned>& ipoints,
                    const Dim<3>::Vector& normal) {
  typedef Dim<3>::Vector Vector;

  // Prerequisites.
  const unsigned npts = ipoints.size();
  unsigned i, j, k;
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(ipoints.size() > 2);
    const double normmag = normal.magnitude();
    Vector normi;
    for (i = 1; i != npts; ++i) {
      j = (i + 1) % npts;
      k = (i + 2) % npts;
      normi = (vertices[ipoints[j]] - vertices[ipoints[i]]).cross(vertices[ipoints[k]] - vertices[ipoints[i]]);
      REQUIRE(fuzzyEqual(abs(normi.dot(normal)), normmag*normi.magnitude(), 1.0e-10));
    }
    REQUIRE(fuzzyEqual(pointPlaneDistance(p, vertices[ipoints[0]], normal.unitVector()), 0.0, 1.0e-10));
  }
  END_CONTRACT_SCOPE;

  bool result = false;

  const double px = p.x();
  const double py = p.y();
  const double pz = p.z();

  double fxmin =  1e50, fymin =  1e50, fzmin =  1e50;
  double fxmax = -1e50, fymax = -1e50, fzmax = -1e50;
  for (i = 0; i != npts; ++i) {
    fxmin = min(fxmin, vertices[ipoints[i]].x());
    fymin = min(fymin, vertices[ipoints[i]].y());
    fzmin = min(fzmin, vertices[ipoints[i]].z());
    fxmax = max(fxmax, vertices[ipoints[i]].x());
    fymax = max(fymax, vertices[ipoints[i]].y());
    fzmax = max(fzmax, vertices[ipoints[i]].z());
  }
  if (px >= fxmin and px <= fxmax and
      py >= fymin and py <= fymax and
      pz >= fzmin and pz <= fzmax) {

    // Check if the point is on the boundary.
    for (i = 0; i != npts; ++i) {
      j = (i + 1) % npts;
      if ((closestPointOnSegment(p, vertices[ipoints[i]], vertices[ipoints[j]]) - p).magnitude2() < 1.0e-10) return true;
    }

    if (abs(normal.x()) >= abs(normal.y()) and abs(normal.x()) >= abs(normal.z())) {

      // x plane -- use (y,z) coordinates.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[ipoints[i]].z() > pz) != (vertices[ipoints[j]].z() > pz)) &&
             (py < (vertices[ipoints[j]].y() - vertices[ipoints[i]].y()) * (pz - vertices[ipoints[i]].z()) / (vertices[ipoints[j]].z() - vertices[ipoints[i]].z()) + vertices[ipoints[i]].y()) )
          result = not result;
      }

    } else if (abs(normal.y()) >= abs(normal.x()) and abs(normal.y()) >= abs(normal.z())) {

      // y plane -- use (x,z) coordinates.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[ipoints[i]].z() > pz) != (vertices[ipoints[j]].z() > pz)) &&
             (px < (vertices[ipoints[j]].x() - vertices[ipoints[i]].x()) * (pz - vertices[ipoints[i]].z()) / (vertices[ipoints[j]].z() - vertices[ipoints[i]].z()) + vertices[ipoints[i]].x()) )
          result = not result;
      }

    } else {
      CHECK(abs(normal.z()) >= abs(normal.x()) and abs(normal.z()) >= abs(normal.y()));

      // z plane -- use (x,y) coordinate.
      for (i = 0, j = npts - 1; i < npts; j = i++) {
        if ( ((vertices[ipoints[i]].y() > py) != (vertices[ipoints[j]].y() > py)) &&
             (px < (vertices[ipoints[j]].x() - vertices[ipoints[i]].x()) * (py - vertices[ipoints[i]].y()) / (vertices[ipoints[j]].y() - vertices[ipoints[i]].y()) + vertices[ipoints[i]].x()) )
          result = not result;
      }

    }
  }

  // That's it.
  return result;
}

}
