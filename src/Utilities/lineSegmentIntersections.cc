//---------------------------------*-C++-*------------------------------------//
// Copyright 2009 LLNS
// All rights reserved.
//
// Created on: Wed Aug 19 10:51:19 PDT 2009
// Created by: J. Michael Owen
// Modified by: 
//
// Methods for computing intersections of line segments with other geometric
// types.
// These are largely adapted from code in "Computational Geometry in C" by
// Joseph O'Rourke.
//----------------------------------------------------------------------------//
#include "lineSegmentIntersections.hh"
#include "pointInPolygon.hh"
#include "spheralWildMagicConverters.hh"
#include "Wm5Segment2.h"
#include "Wm5Segment3.h"
#include "Wm5DistSegment2Segment2.h"
#include "Wm5DistSegment3Segment3.h"

using namespace std;

namespace Spheral {

//------------------------------------------------------------------------------
// A utility to make us scale independent by determining how to scale input
// points.
//------------------------------------------------------------------------------
template<typename Vector>
inline
double
helpfulScaleFactor(const Vector& a,
                   const Vector& b,
                   const Vector& c) {
  return std::max(1e-10, std::max(a.maxAbsElement(), std::max(b.maxAbsElement(), c.maxAbsElement())));
}

template<typename Vector>
inline
double
helpfulScaleFactor(const Vector& a,
                   const Vector& b,
                   const Vector& c,
                   const Vector& d) {
  return std::max(helpfulScaleFactor(a, b, c), d.maxAbsElement());
}

//------------------------------------------------------------------------------
// Test for intersection between two parallel segments.
//------------------------------------------------------------------------------
template<typename Vector>
char
parallelSegmentIntersection(const Vector& a0,
                            const Vector& a1,
                            const Vector& b0,
                            const Vector& b1,
                            Vector& result1,
                            Vector& result2,
                            const double tol) {
  REQUIRE(fuzzyEqual(abs((a1 - a0).dot(b1 - b0)), (a1 - a0).magnitude() * (b1 - b0).magnitude(), tol));
  if (not collinear(a0, a1, b0, tol)) {
    result1 = Vector::zero;
    result2 = Vector::zero;
    return '0';
  }
  const bool a0test = between(b0, b1, a0, tol);
  const bool a1test = between(b0, b1, a1, tol);
  const bool b0test = between(a0, a1, b0, tol);
  const bool b1test = between(a0, a1, b1, tol);
  const unsigned num = ((a0test ? 1 : 0) +
                        (a1test ? 1 : 0) +
                        (b0test ? 1 : 0) +
                        (b1test ? 1 : 0));
  CHECK(num <= 4);

  // No intersection.
  if (num == 0) {
    result1 = Vector::zero;
    result2 = Vector::zero;
    return '0';
  }

  // One intersection.
  if (num == 1) {
    if (a0test) {
      result1 = a0;
      result2 = a0;
      return 'e';
    }
    if (a1test) {
      result1 = a1;
      result2 = a1;
      return 'e';
    }
    if (b0test) {
      result1 = b0;
      result2 = b0;
      return 'e';
    }
    if (b1test) {
      result1 = b1;
      result2 = b1;
      return 'e';
    }
    CHECK(false);
  }

  // Two or more intersections.
  if (a0test and a1test) {
    result1 = a0;
    result2 = a1;
    return 'e';
  }
  if (b0test and b1test) {
    result1 = b0;
    result2 = b1;
    return 'e';
  }
  if (a0test and b0test) {
    result1 = a0;
    result2 = a1;
    return 'e';
  }
  if (a0test and b1test) {
    result1 = a0;
    result2 = b1;
    return 'e';
  }
  if (a1test and b0test) {
    result1 = a1;
    result2 = b0;
    return 'e';
  }
  if (a1test and b1test) {
    result1 = a1;
    result2 = b1;
    return 'e';
  }
  CHECK(false);
}

//------------------------------------------------------------------------------
// Intersection of two line segments.
// The line segments are characterized by their endpoints: a_seg = (a0, a1)
//                                                         b_seg = (b0, b1)
// Return values are a pair<Vector, char>
//        The Vector is the intersection point (if any)
//        The char is a code characterizing the intersection:
//              "e" -> The segments colinearly overlap (edge)
//              "v" -> The endpoint of one segment lies on the other
//              "1" -> The segments intersect properly
//              "0" -> The segments do not intersect
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
char
segmentSegmentIntersection(const Dim<1>::Vector& a0,
                           const Dim<1>::Vector& a1,
                           const Dim<1>::Vector& b0,
                           const Dim<1>::Vector& b1,
                           Dim<1>::Vector& result1,
                           Dim<1>::Vector& result2,
                           const double tol) {
  typedef Dim<1>::Vector Vector;
  const double reltol = tol*helpfulScaleFactor(a0, a1, b0, b1);

  vector<Vector> results;
  if (between(a0, a1, b0, reltol)) results.push_back(b0);
  if (between(a0, a1, b1, reltol)) results.push_back(b1);
  if (between(b0, b1, a0, reltol)) results.push_back(a0);
  if (between(b0, b1, a1, reltol)) results.push_back(a1);
  CHECK(results.size() <= 2);

  if (results.size() == 0) {
    return '0';
  } else if (results.size() == 1) {
    result1 = results[0];
    result2 = result1;
    return 'e';
  } else {
    result1 = results[0];
    result2 = results[1];
    return 'e';
  }
}    

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
char
segmentSegmentIntersection(const Dim<2>::Vector& a0,
                           const Dim<2>::Vector& a1,
                           const Dim<2>::Vector& b0,
                           const Dim<2>::Vector& b1,
                           Dim<2>::Vector& result1,
                           Dim<2>::Vector& result2,
                           const double tol) {
  typedef Dim<2>::Vector Vector;
  const double reltol = tol*helpfulScaleFactor(a0, a1, b0, b1);

  // Prepare the results.
  char code = 'z';

  const double denom = (a0.x()*(b1.y() - b0.y()) +
                        a1.x()*(b0.y() - b1.y()) +
                        b1.x()*(a1.y() - a0.y()) +
                        b0.x()*(a0.y() - a1.y()));

  // If the denominator is 0, the segments are parallel.
  if (fuzzyEqual(denom, 0.0, reltol)) return parallelSegmentIntersection(a0, a1, b0, b1, result1, result2, reltol);

  const double num1 = (a0.x()*(b1.y() - b0.y()) +
                       b0.x()*(a0.y() - b1.y()) +
                       b1.x()*(b0.y() - a0.y()));
  const double num2 = -(a0.x()*(b0.y() - a1.y()) +
                        a1.x()*(a0.y() - b0.y()) +
                        b0.x()*(a1.y() - a0.y()));
  const double s = num1/denom;
  const double t = num2/denom;

  // Determine the type of intersection.
  if (between(b0, b1, a0, reltol) or
      between(b0, b1, a1, reltol) or
      between(a0, a1, b0, reltol) or 
      between(a0, a1, b1, reltol)) {

    // Endpoint of one segment on the other.
    code = 'v';

  } else if (s >= 0.0 and s <= 1.0 and t >= 0.0 and t <= 1.0) {

    // Proper intersection.
    code = '1';

  } else if (s < 0.0 or s > 1.0 or t < 0.0 or t > 1.0) {

    // No intersection.
    code = '0';

  }

  result1 = a0 + s*(a1 - a0);
  result2 = result1;
  ENSURE(code != 'z');
  return code;
}

//------------------------------------------------------------------------------
// 3D.
//------------------------------------------------------------------------------
char
segmentSegmentIntersection(const Dim<3>::Vector& a0,
                           const Dim<3>::Vector& a1,
                           const Dim<3>::Vector& b0,
                           const Dim<3>::Vector& b1,
                           Dim<3>::Vector& result1,
                           Dim<3>::Vector& result2,
                           const double tol) {
  typedef Dim<3>::Vector Vector;
  const double reltol = tol*helpfulScaleFactor(a0, a1, b0, b1);

  // Are the segments parallel?
  const Vector na = (a1 - a0).unitVector();
  const Vector nb = (b1 - b0).unitVector();
  if (fuzzyEqual(std::abs(na.dot(nb)), 1.0, tol)) return parallelSegmentIntersection(a0, a1, b0, b1, result1, result2, reltol);

  // Wild magic has accuracy problems when the endpoint of one segment lies on 
  // another, so we screen for that case.
  if (between(a0, a1, b0, reltol)) { result1 = b0; result2 = b0; return 'v'; }
  if (between(a0, a1, b1, reltol)) { result1 = b1; result2 = b1; return 'v'; }
  if (between(b0, b1, a0, reltol)) { result1 = a0; result2 = a0; return 'v'; }
  if (between(b0, b1, a1, reltol)) { result1 = a1; result2 = a1; return 'v'; }

  // Let WildMagic compute the general case.
  const Wm5::Segment3<double> aseg(convertVectorToWMVector<Dim<3> >(a0),
                                   convertVectorToWMVector<Dim<3> >(a1));
  const Wm5::Segment3<double> bseg(convertVectorToWMVector<Dim<3> >(b0),
                                   convertVectorToWMVector<Dim<3> >(b1));
  Wm5::DistSegment3Segment3<double> dist(aseg, bseg);
  if (fuzzyEqual(dist.Get(), 0.0, reltol)) {
    result1 = a0 + (a1 - a0)*dist.GetSegment0Parameter();
    result2 = result1;
    return '1';
  }
}

//------------------------------------------------------------------------------
// Intersection of a line segment with a plane.
// The line segment is characterized by it's endpoints:     seg = (s0, s1)
// The plane is characterized by three points in the plane: plane = (p0, p1, p2)
//
// Return values are a pair<Vector, char>
//        The Vector is the intersection point (if any)
//        The char is a code characterizing the intersection:
//              "p" -> The segment lies in the plane (plane)
//              "d" -> The p points do not define a unique plane (degenerate)
//              "1" -> The segment intersects the plane properly
//              "0" -> The segment does not intersect the plane
//------------------------------------------------------------------------------
char
segmentPlaneIntersection(const Dim<3>::Vector& s0,
                         const Dim<3>::Vector& s1,
                         const Dim<3>::Vector& p,
                         const Dim<3>::Vector& n,
                         Dim<3>::Vector& result,
                         const double tol) {
  typedef Dim<3>::Vector Vector;

  // Determine the equation of the plane.
  const double reltol = tol*helpfulScaleFactor(s0, s1, p);
  const double D = n.dot(p);

  // If the normal is zero, then the plane points are degenerate and we're done.
  if (fuzzyEqual(n.magnitude(), 0.0, reltol)) {
    result.Zero();
    return 'd';
  }

  // Is the line segment in the plane?
  if (fuzzyEqual(n.dot(s0) - D, 0.0, reltol) and fuzzyEqual(n.dot(s1) - D, 0.0, reltol)) {
    result = s0;
    return 'p';
  }

  // Is the line segment parallel to the plane?
  const Vector s10 = s1 - s0;
  const double den = n.dot(s10);
  if (fuzzyEqual(den, 0.0, reltol)) {
    result.Zero();
    return '0';
  }

  // Compute the intersection.  If q is outside the range [0,1], then the line
  // segment does not intersect the plane.
  const double q = (D - n.dot(s0))/den;
  if (fuzzyGreaterThanOrEqual(q, 0.0, tol) and
      fuzzyLessThanOrEqual(q, 1.0, tol)) {
    result = s0 + q*s10;
    return '1';
  } else {
    result.Zero();
    return '0';
  }
}

//------------------------------------------------------------------------------
// Intersection of a line segment with a plane.
// The line segment is characterized by it's endpoints:     seg = (s0, s1)
// The plane is characterized by three points in the plane: plane = (p0, p1, p2)
//
// Return values are a pair<Vector, char>
//        The Vector is the intersection point (if any)
//        The char is a code characterizing the intersection:
//              "p" -> The segment lies in the plane (plane)
//              "d" -> The p points do not define a unique plane (degenerate)
//              "1" -> The segment intersects the plane properly
//              "0" -> The segment does not intersect the plane
//------------------------------------------------------------------------------
char
segmentPlaneIntersection(const Dim<3>::Vector& s0,
                         const Dim<3>::Vector& s1,
                         const Dim<3>::Vector& p0,
                         const Dim<3>::Vector& p1,
                         const Dim<3>::Vector& p2,
                         Dim<3>::Vector& result,
                         const double tol) {
  typedef Dim<3>::Vector Vector;

  // Determine the equation of the plane.
  const Vector p10 = p1 - p0;
  const Vector p20 = p2 - p0;
  const Vector n = p20.cross(p10);

  return segmentPlaneIntersection(s0, s1, p0, n, result, tol);
}

//------------------------------------------------------------------------------
// Intersection of a line segment with a polygonal section of a plane.
// The line segment is characterized by it's endpoints:     seg = (s0, s1)
// The polygonal section of the plane is specified by
// a series of points:                                    plane = pverts[0], pverts[1], ...
// Note there must be at least 3 non-collinear points, and they
// must be passed in in order to draw the polygon.
//
// Return values are a pair<Vector, char>
//        The Vector is the intersection point (if any)
//        The char is a code characterizing the intersection:
//              "p" -> The segment lies in the plane (plane) and the segment hits it.
//              "d" -> The pverts do not define a unique plane (degenerate)
//              "1" -> The segment intersects the plane properly and is in the polygon
//              "0" -> The segment intersection does not lie in the polygon
//------------------------------------------------------------------------------
char
segmentPlanarSectionIntersection(const Dim<3>::Vector& s0,
                                 const Dim<3>::Vector& s1,
                                 const vector<Dim<3>::Vector>& pverts,
                                 Dim<3>::Vector& result,
                                 const double tol) {
  typedef Dim<3>::Vector Vector;
  REQUIRE2(pverts.size() > 2, "Need at least three points to specify a plane.");

  // Find three non-degenerate points to define our plane.
  unsigned ip1 = 0;
  while (ip1 < pverts.size() and (pverts[ip1] - pverts[0]).magnitude2() <= tol) ++ip1;
  if (ip1 == pverts.size()) return 'd';
  unsigned ip2 = ip1 + 1;
  while (ip2 < pverts.size() and collinear(pverts[0], pverts[ip1], pverts[ip2])) ++ip2;
  if (ip2 == pverts.size()) return 'd';

  // Determine if the point intersects the plane.
  char code = segmentPlaneIntersection(s0, s1, 
                                       pverts[0], pverts[ip1], pverts[ip2],
                                       result, tol);
  CHECK(code != 'd');

  // If the segment doesn't intersect the plane, we're done.
  char code2;
  if (code == '0') {
    return code;

  } else if (code == '1') {
    // The segment intersects the plane.  Check if that intersection is in the polygon.
    const Vector normal = (pverts[ip1] - pverts[0]).cross(pverts[ip2] - pverts[0]);
    if (pointInPolygon(result, pverts, normal)) { 
      return '1';
    } else {
      return '0';
    }

  } else {
    CHECK(code == 'p');
    // The segment lies in the plane of the polygon.  This case is not unique,
    // so check if either endpoint is in the polygon, or if the segment intersects
    // any of the polygon's edges.
    const Vector normal = (pverts[ip1] - pverts[0]).cross(pverts[ip2] - pverts[0]);
    if (pointInPolygon(s0, pverts, normal)) {
      result = s0;
      return 'p';
    } else if (pointInPolygon(s1, pverts, normal)) {
      result = s1;
      return 'p';
    } else {
      unsigned i = 0, j;
      char code2;
      while (i < pverts.size()) {
        j = i % pverts.size();
        code2 = segmentSegmentIntersection(s0, s1, pverts[i], pverts[j], result, result, tol);
        if (code2 != '0') return 'p';
        ++i;
      }
    }
  }

  return '0';
}

//------------------------------------------------------------------------------
// Compute the minimum distance between two line segments. (a0,a1) -> (b0,b1).
// In 2-D and 3-D we simply farm this work out to David Eberly's code from 
// www.geometrictools.com
//------------------------------------------------------------------------------
// 1-D.
double
segmentSegmentDistance(const Dim<1>::Vector& a0,
                       const Dim<1>::Vector& a1,
                       const Dim<1>::Vector& b0,
                       const Dim<1>::Vector& b1) {

  const double amin = std::min(a0.x(), a1.x());
  const double amax = std::max(a0.x(), a1.x());
  const double bmin = std::min(b0.x(), b1.x());
  const double bmax = std::max(b0.x(), b1.x());

  // Do the segments overlap?
  if ((amin <= bmin and bmin <= amax) or
      (amin <= bmax and bmax <= amax) or
      (bmin <= amin and amin <= bmax) or
      (bmin <= amax and amax <= bmax)) return 0.0;

  // Otherwise it's the minimum distance between the endpoints.
  return std::min(std::abs(bmin - amax), std::abs(amin - bmax));
}

// 2-D.
double
segmentSegmentDistance(const Dim<2>::Vector& a0,
                       const Dim<2>::Vector& a1,
                       const Dim<2>::Vector& b0,
                       const Dim<2>::Vector& b1) {
  // Wild magic has accuracy problems when the endpoint of one segment lies on 
  // another, so we screen for that case.
  const double reltol = 1.0e-8*helpfulScaleFactor(a0, a1, b0, b1);
  if (between(a0, a1, b0, reltol) or between(a0, a1, b1, reltol) or
      between(b0, b1, a0, reltol) or between(b0, b1, a1, reltol)) return 0.0;

  const Wm5::Segment2<double> aseg(convertVectorToWMVector<Dim<2> >(a0),
                                   convertVectorToWMVector<Dim<2> >(a1));
  const Wm5::Segment2<double> bseg(convertVectorToWMVector<Dim<2> >(b0),
                                   convertVectorToWMVector<Dim<2> >(b1));
  return Wm5::DistSegment2Segment2<double>(aseg, bseg).Get();
}

// 3-D.
double
segmentSegmentDistance(const Dim<3>::Vector& a0,
                       const Dim<3>::Vector& a1,
                       const Dim<3>::Vector& b0,
                       const Dim<3>::Vector& b1) {
  // Wild magic has accuracy problems when the endpoint of one segment lies on 
  // another, so we screen for that case.
  const double reltol = 1.0e-8*helpfulScaleFactor(a0, a1, b0, b1);
  if (between(a0, a1, b0, reltol) or between(a0, a1, b1, reltol) or
      between(b0, b1, a0, reltol) or between(b0, b1, a1, reltol)) return 0.0;

  const Wm5::Segment3<double> aseg(convertVectorToWMVector<Dim<3> >(a0),
                                   convertVectorToWMVector<Dim<3> >(a1));
  const Wm5::Segment3<double> bseg(convertVectorToWMVector<Dim<3> >(b0),
                                   convertVectorToWMVector<Dim<3> >(b1));
  return Wm5::DistSegment3Segment3<double>(aseg, bseg).Get();
}

//------------------------------------------------------------------------------
// Check if a line segment intersects a polygonal planar section in 3D.
// The line segment is characterized by its endpoints: a_seg = (a0, a1)
// The planar section is characterized the a set of coplanar vertices: vertices.
//
// The advantage of this method is that it should be faster if you don't need
// the actual intersection.
// 
// The second version allows the planar vertices to be represented as a subset
// of the set passed in.
//------------------------------------------------------------------------------
bool
segmentPlaneIntersection(const Dim<3>::Vector& a0,
                         const Dim<3>::Vector& a1,
                         const std::vector<Dim<3>::Vector>& vertices,
                         const Dim<3>::Vector& normal,
                         const double tol) {
  REQUIRE(vertices.size() > 2);
  typedef Dim<3>::Vector Vector;
  Vector p;
  const char code = segmentPlaneIntersection(a0, a1, vertices[0], normal, p, tol);
  if (code == '0') return false;
  if (code == '1') return pointInPolygon(p, vertices, normal);

  // The line segment is in the plane -- this case is degenerate.
  CHECK(code == 'p');
  if (pointInPolygon(a0, vertices, normal) or pointInPolygon(a1, vertices, normal)) return true;
  unsigned i, j, n = vertices.size();
  for (i = 0; i != n; ++i) {
    j = (i + 1) % n;
    if (segmentSegmentIntersection(a0, a1, vertices[i], vertices[j], tol)) return true;
  }
  return false;
}

bool
segmentPlaneIntersection(const Dim<3>::Vector& a0,
                         const Dim<3>::Vector& a1,
                         const std::vector<Dim<3>::Vector>& vertices,
                         const std::vector<unsigned>& ipoints,
                         const Dim<3>::Vector& normal,
                         const double tol) {
  REQUIRE(vertices.size() > 2);
  REQUIRE(ipoints.size() > 2);
  REQUIRE(vertices.size() >= ipoints.size());
  typedef Dim<3>::Vector Vector;
  Vector p;
  const char code = segmentPlaneIntersection(a0, a1, vertices[ipoints[0]], normal, p, tol);
  if (code == '0') return false;
  if (code == '1') return pointInPolygon(p, vertices, ipoints, normal);

  // The line segment is in the plane -- this case is degenerate.
  CHECK(code == 'p');
  if (pointInPolygon(a0, vertices, ipoints, normal) or pointInPolygon(a1, vertices, ipoints, normal)) return true;
  unsigned i, j, n = ipoints.size();
  for (i = 0; i != n; ++i) {
    j = (i + 1) % n;
    if (segmentSegmentIntersection(a0, a1, vertices[ipoints[i]], vertices[ipoints[j]], tol)) return true;
  }
  return false;
}


}
