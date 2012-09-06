//---------------------------------Spheral++----------------------------------//
// GeomPolygon -- Geometric polygon class.
//
// A 2-D structure representing a polygon as a collection of GeomFacets.
//
// Created by JMO, Thu Jan 28 11:03:27 PST 2010
//----------------------------------------------------------------------------//
#include <algorithm>
#include <numeric>
#include <map>
#include <limits>

#include "boost/foreach.hpp"

#include "polytope/polytope.hh"
#include "polytope/convexHull_2d.hh"

#include "GeomPolygon.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/spheralWildMagicConverters.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/pointInPolygon.hh"

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
GeomPolygon::
GeomPolygon():
  mVertices(),
  mFacets(),
  mXmin(),
  mXmax(),
  mConvex(true) {
  if (mDevnull == NULL) mDevnull = fopen("/dev/null", "w");
}

//------------------------------------------------------------------------------
// Construct as a convex hull for the given point set.
//------------------------------------------------------------------------------
GeomPolygon::
GeomPolygon(const vector<GeomPolygon::Vector>& points):
  mVertices(),
  mFacets(),
  mXmin(),
  mXmax(),
  mConvex(true) {
  if (mDevnull == NULL) mDevnull = fopen("/dev/null", "w");

  if (points.size() > 0) {

    REQUIRE(points.size() > 2);

    // Find the appropriate renormalization so that we can do the convex hull
    // in a unit box.
    Vector xmin, xmax;
    boundingBox(points, xmin, xmax);
    const double fscale = (xmax - xmin).maxElement();
    CHECK(fscale > 0.0);

    // Copy the point coordinates to a polytope point array.
    vector<double> points_polytope;
    points_polytope.reserve(2 * points.size());
    BOOST_FOREACH(Vector vec, points) {
      points_polytope.push_back((vec.x() - xmin.x())/fscale);
      points_polytope.push_back((vec.y() - xmin.y())/fscale);
    }
    CHECK(points_polytope.size() == 2*points.size());

    // Call the polytope method for computing the convex hull.
    vector<double> low(2, 0.0);
    polytope::PLC<2, double> plc = polytope::convexHull_2d(points_polytope, &(*low.begin()), 1.0e-10);
    const unsigned numVertices = plc.facets.size();
    CHECK(numVertices >= 3);

    // Extract the hull information back to our local convention.  We use the fact that
    // polytope's convex hull method sorts the vertices in counter-clockwise here.
    // Start with the vertices.
    mVertices.reserve(numVertices);
    int i, j;
    for (j = 0; j != numVertices; ++j) {
      CHECK(plc.facets[j].size() == 2);
      i = plc.facets[j][0];
      CHECK(i >= 0 and i < points.size());
      mVertices.push_back(points[i]);
    }

    // Now the facets.
    mFacets.reserve(numVertices);
    for (i = 0; i != numVertices; ++i) {
      j = (i + 1) % numVertices;
      mFacets.push_back(Facet(mVertices, i, j));
    }

    // Fill in our bounding box.
    setBoundingBox();

    // Post-conditions.
    BEGIN_CONTRACT_SCOPE;
    {
      // Ensure the facet node ordering is correct.
      CounterClockwiseComparator<Vector, vector<Vector> > nodeComparator(mVertices, mVertices[0]);
      BOOST_FOREACH(const Facet& facet, mFacets) ENSURE(nodeComparator(facet.point1(), facet.point2()) >= 0);

      // All normals should be outward facing.
      Vector centroid, vec;
      BOOST_FOREACH(vec, mVertices) centroid += vec;
      centroid /= mVertices.size();
      BOOST_FOREACH(const Facet& facet, mFacets) ENSURE2((0.5*(facet.point1() + facet.point2()) - centroid).dot(facet.normal()) >= 0.0,
                                                         facet.point1() << " " << facet.point2() << " : "
                                                         << (0.5*(facet.point1() + facet.point2()) - centroid) << " "
                                                         << facet.normal() << " : "
                                                         << (0.5*(facet.point1() + facet.point2()) - centroid).dot(facet.normal()));

      // Ensure the vertices are listed in counter-clockwise order.
      for (unsigned i = 0; i != mVertices.size(); ++i) {
        const unsigned j = (i + 1) % mVertices.size();
        ENSURE(nodeComparator(i, j) >= 0);
      }

      // We had better be convex if built from a convex hull.
      ENSURE(this->convex());

      // Ensure the seed points are contained.
      // Suspending this check for now as floating point accuracy occasionally misfires
      // this check.
      //      BOOST_FOREACH(vec, points) ENSURE(this->convexContains(vec));
    }
    END_CONTRACT_SCOPE;
  }
}

//------------------------------------------------------------------------------
// Construct given the positions and facet indices.
// We assume here that the nodes for each facet are arranged correctly to
// create outward pointing normals.
//------------------------------------------------------------------------------
GeomPolygon::
GeomPolygon(const vector<GeomPolygon::Vector>& points,
            const vector<vector<unsigned> >& facetIndices):
  mVertices(points),
  mFacets(),
  mConvex(false) {

  // Construct the facets.
  Vector centroid;
  mFacets.reserve(facetIndices.size());
  BOOST_FOREACH(vector<unsigned> indices, facetIndices) {
    VERIFY2(indices.size() == 2, "Need two points per facet : " << indices.size());
    VERIFY2(*max_element(indices.begin(), indices.end()) < points.size(),
            "Bad vertex index for facet.");
    mFacets.push_back(Facet(mVertices, indices[0], indices[1]));
  }
  CHECK(mFacets.size() == facetIndices.size());

  // Fill in our bounding box.
  setBoundingBox();

  // Check if we're convex.
  mConvex = this->convex();
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
GeomPolygon::
GeomPolygon(const GeomPolygon& rhs):
  mVertices(rhs.mVertices),
  mFacets(),
  mXmin(rhs.mXmin),
  mXmax(rhs.mXmax),
  mConvex(rhs.mConvex) {
  mFacets.reserve(rhs.mFacets.size());
  BOOST_FOREACH(const Facet& facet, rhs.mFacets) {
    mFacets.push_back(Facet(mVertices,
                            facet.ipoint1(),
                            facet.ipoint2()));
  }
  ENSURE(mFacets.size() == rhs.mFacets.size());
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
GeomPolygon&
GeomPolygon::
operator=(const GeomPolygon& rhs) {
  if (this != &rhs) {
    mVertices = rhs.mVertices;
    mFacets = vector<Facet>();
    mFacets.reserve(rhs.mFacets.size());
    BOOST_FOREACH(const Facet& facet, rhs.mFacets) {
         mFacets.push_back(Facet(mVertices,
                                 facet.ipoint1(),
                                 facet.ipoint2()));
    }
    mXmin = rhs.mXmin;
    mXmax = rhs.mXmax;
    mConvex = rhs.mConvex;
  }
  ENSURE(mFacets.size() == rhs.mFacets.size());
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
GeomPolygon::
~GeomPolygon() {
}

//------------------------------------------------------------------------------
// Test if the given point is on the interior of the polygon.
// Generic test.
//------------------------------------------------------------------------------
bool
GeomPolygon::
contains(const GeomPolygon::Vector& point,
         const bool countBoundary,
         const double tol) const {
  if (mConvex) {
    cerr << "GeomPolygon : executing convex test." << endl;
    return this->convexContains(point, countBoundary, tol);
  } else {
    cerr << "GeomPolygon : executing non-convex test." << endl;
    return pointInPolygon(point, *this, countBoundary, tol);
  }
}

//------------------------------------------------------------------------------
// Test if the given point is on the interior of the polygon.
// This method only works for convex polygons!
//------------------------------------------------------------------------------
bool
GeomPolygon::
convexContains(const GeomPolygon::Vector& point,
               const bool countBoundary,
               const double tol) const {
  if (not testPointInBox(point, mXmin, mXmax, tol)) return false;
  vector<Facet>::const_iterator facetItr = mFacets.begin();
  bool result = true;
  if (countBoundary) {
    while (facetItr != mFacets.end() and result) {
      result = (facetItr->compare(point, tol) <= 0);
      if (not result) {
        cerr << "Failing facet test : " << *facetItr << " " << point << endl;
      }
      ++facetItr;
    }
  } else {
    while (facetItr != mFacets.end() and result) {
      result = (facetItr->compare(point, tol) < 0);
      if (not result) {
        cerr << "Failing facet test : " << *facetItr << " " << point << endl;
      }
      ++facetItr;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Test if we intersect the given polygon.
//------------------------------------------------------------------------------
bool
GeomPolygon::
intersect(const GeomPolygon& rhs) const {
  if (not testBoxIntersection(mXmin, mXmax, rhs.mXmin, rhs.mXmax)) return false;
  Vector vec;
  BOOST_FOREACH(vec, mVertices) {
    if (rhs.contains(vec)) return true;
  }
  BOOST_FOREACH(vec, rhs.mVertices) {
    if (this->contains(vec)) return true;
  }
  unsigned i0, j0, i1, j1;
  const unsigned n0 = mVertices.size();
  const unsigned n1 = rhs.mVertices.size();
  for (i0 = 0; i0 != n0; ++i0) {
    j0 = (i0 + 1) % n0;
    for (i1 = 0; i1 != n1; ++i1) {
      j1 = (i1 + 1) % n1;
      if (segmentSegmentIntersection(mVertices[i0], mVertices[j0],
                                     rhs.mVertices[i1], rhs.mVertices[j1])) return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------
// Test if we intersect the given polygon using the knowledge that both polygons
// are convex.
// We use the method of separating axes here.
//------------------------------------------------------------------------------
bool
GeomPolygon::
convexIntersect(const GeomPolygon& rhs) const {
  REQUIRE(this->convex());
  if (not testBoxIntersection(mXmin, mXmax, rhs.mXmin, rhs.mXmax)) return false;
  
  // Check if we can exclude rhs from us.
  bool outside = true;
  {
    std::vector<Facet>::const_iterator facetItr = mFacets.begin();
    while (outside and facetItr != mFacets.end()) {
      outside = (facetItr->compare(rhs.mVertices) == 1);
      ++facetItr;
    }
    if (outside) return false;
  }

  // Check if we can exclude us from rhs.
  outside = true;
  {
    std::vector<Facet>::const_iterator facetItr = rhs.mFacets.begin();
    while (outside and facetItr != rhs.mFacets.end()) {
      outside = (facetItr->compare(mVertices) == 1);
      ++facetItr;
    }
    if (outside) return false;
  }

  // We can't exclude anybody, so must intersect!
  return true;
}

//------------------------------------------------------------------------------
// Test if we intersect a box.
//------------------------------------------------------------------------------
bool
GeomPolygon::
intersect(const GeomPolygon::Box& rhs) const {
  BOOST_FOREACH(Vector vec, mVertices) {
    if (testPointInBox(vec, rhs)) return true;
  }

  typedef Wm5::Vector2<double> WMVector;
  vector<WMVector> WMvertices(4);
  rhs.ComputeVertices(&WMvertices.front());
  if (this->contains(convertWMVectorToVector<Dim<2> >(WMvertices[0]))) return true;
  if (this->contains(convertWMVectorToVector<Dim<2> >(WMvertices[1]))) return true;
  if (this->contains(convertWMVectorToVector<Dim<2> >(WMvertices[2]))) return true;
  if (this->contains(convertWMVectorToVector<Dim<2> >(WMvertices[3]))) return true;

  return false;
}

//------------------------------------------------------------------------------
// Compute the intersections with a line segment.
//------------------------------------------------------------------------------
vector<GeomPolygon::Vector>
GeomPolygon::
intersect(const GeomPolygon::Vector& x0, const GeomPolygon::Vector& x1) const {
  vector<Vector> result;
  Vector inter1, inter2;
  for (size_t inode = 0; inode != mVertices.size(); ++inode) {
    const Vector& e0 = mVertices[inode];
    const Vector& e1 = mVertices[inode % mVertices.size()];
    const char code = segmentSegmentIntersection(x0, x1, e0, e1, inter1, inter2);

    if (code == '1' or code == 'v') {
      // Proper intersection.
      result.push_back(inter1);

    } else if (code == 'e') {
      // The segment is colinear with and overlaps the edge.  In this case we 
      // return the end-points of the section that is on both segments.
      result.push_back(inter1);
      result.push_back(inter2);
      // if (between(x0, x1, e0)) result.push_back(e0);
      // if (between(x0, x1, e1)) result.push_back(e1);
      // if (between(e0, e1, x0)) result.push_back(x0);
      // if (between(e0, e1, x1)) result.push_back(x1);
    }
  }

  // It's possible that we may have duplicates in the intersection set. 
  // Make it unique!
  sort(result.begin(), result.end());
  result.erase(unique(result.begin(), result.end()), result.end());

  // That's it.
  ENSURE(result.size() <= 2);
  return result;
}

//------------------------------------------------------------------------------
// Compute the centroid.
//------------------------------------------------------------------------------
GeomPolygon::Vector
GeomPolygon::
centroid() const {
  const unsigned n = mVertices.size();
  Vector result;
  if (n == 0) return result;
  CHECK(n >= 3);
  const Vector c0 = std::accumulate(mVertices.begin(), mVertices.end(), Vector())/n;

  // Specialize for triangles, which are easy!
  if (n == 3) return c0;

  unsigned i, j;
  double area, areasum = 0.0;
  for (i = 0; i != n; ++i) {
    j = (i + 1) % n;
    area = (mVertices[i] - c0).cross(mVertices[j] - c0).z(); // This is off by a factor of 2 but will cancel.
    areasum += area;
    result += area * (c0 + mVertices[i] + mVertices[j]);
  }
  CHECK(areasum > 0.0);
  return result/(3.0 * areasum);
}

//------------------------------------------------------------------------------
// Return the edges of the polygon as a set of integer pairs for the vertices.
// This is simplified because we know the vertices are already sorted 
// counter-clockwise.
//------------------------------------------------------------------------------
vector<pair<unsigned, unsigned> >
GeomPolygon::
edges() const {
  vector<pair<unsigned, unsigned> > result;
  for (unsigned i = 0; i != mVertices.size(); ++i) {
    result.push_back(make_pair(i, (i + 1) % mVertices.size()));
  }
  return result;
}

//------------------------------------------------------------------------------
// Spit out an encoding of the facets as ordered vertex indices.
//------------------------------------------------------------------------------
vector<vector<unsigned> >
GeomPolygon::
facetVertices() const {
  vector<vector<unsigned> > result;
  vector<unsigned> pts(2);
  if (mVertices.size() > 0) {
    BOOST_FOREACH(const Facet& facet, mFacets) {
      pts[0] = facet.ipoint1();
      pts[1] = facet.ipoint2();
      CHECK(pts[0] < mVertices.size());
      CHECK(pts[1] < mVertices.size());
      result.push_back(pts);
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Reconstruct the internal state given the set of vertices and the enocded 
// facets.
//------------------------------------------------------------------------------
void
GeomPolygon::
reconstruct(const vector<GeomPolygon::Vector>& vertices,
            const vector<vector<unsigned> >& facetVertices) {
  mVertices = vertices;
  mFacets = vector<Facet>();
  mFacets.reserve(facetVertices.size());
  BOOST_FOREACH(const vector<unsigned>& ipts, facetVertices) {
    CHECK2(ipts.size() == 2, "Bad size:  " << ipts.size());
    mFacets.push_back(Facet(mVertices, ipts[0], ipts[1]));
  }
  setBoundingBox();
  mConvex = this->convex();
  ENSURE(mFacets.size() == facetVertices.size());
}

//------------------------------------------------------------------------------
// Compute the volume.
//------------------------------------------------------------------------------
double
GeomPolygon::
volume() const {
  double result = 0.0;
  const Vector c = centroid();
  BOOST_FOREACH(const Facet& facet, mFacets) {
    result += ((facet.point2() - facet.point1()).cross(c - facet.point1())).z();
  }
  ENSURE2(result >= 0.0, result);
  return 0.5*result;
}

//------------------------------------------------------------------------------
// Find the minimum distance to a point.
//------------------------------------------------------------------------------
double
GeomPolygon::
distance(const GeomPolygon::Vector& p) const {
  return (p - this->closestPoint(p)).magnitude();
}

//------------------------------------------------------------------------------
// Find the point in the polygon closest to the given point.
//------------------------------------------------------------------------------
GeomPolygon::Vector
GeomPolygon::
closestPoint(const GeomPolygon::Vector& p) const {
  double r2, minr2 = numeric_limits<double>::max();
  Vector result, thpt;
  BOOST_FOREACH(const Facet& facet, mFacets) {
    thpt = facet.closestPoint(p);
    r2 = (thpt - p).magnitude2();
    if (r2 < minr2) {
      result = thpt;
      minr2 = r2;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
bool
GeomPolygon::
operator==(const GeomPolygon& rhs) const {
  bool result = (mVertices == rhs.mVertices and
                 mFacets.size() == rhs.mFacets.size());
  int i = 0;
  while (result and i != mFacets.size()) {
    result = mFacets[i] == rhs.mFacets[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
bool
GeomPolygon::
operator!=(const GeomPolygon& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// Set the minimum and maximum extents of the polygon.
//------------------------------------------------------------------------------
void
GeomPolygon::
setBoundingBox() {
  boundingBox(mVertices, mXmin, mXmax);
}

//------------------------------------------------------------------------------
// Test if the polygon is convex.
//------------------------------------------------------------------------------
bool
GeomPolygon::
convex(const double tol) const {
  // Do the convex comparison for each vertex.
  bool result = true;
  const double reltol = tol*max(1.0, (mXmax - mXmin).maxAbsElement());
  vector<Vector>::const_iterator vertexItr = mVertices.begin();
  while (vertexItr != mVertices.end() and result) {
    vector<Facet>::const_iterator facetItr = mFacets.begin();
    while (facetItr != mFacets.end() and result) {
      result = (facetItr->compare(*vertexItr, reltol) <= 0);
      if (not result) {
        cerr << "Convex test failing for " << *facetItr << " " << *vertexItr << endl;
      }
      ++facetItr;
    }
    ++vertexItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// ostream operator.
//------------------------------------------------------------------------------
ostream& operator<<(ostream& os, const GeomPolygon& polygon) {
  typedef GeomPolygon::Vector Vector;
  typedef GeomPolygon::Facet Facet;
  const vector<Vector>& vertices = polygon.vertices();
  if (vertices.size() > 0) {
    for (unsigned i = 0; i != vertices.size(); ++i) {
      os << vertices[i].x() << " " << vertices[i].y() << endl;
    }
  }
  return os;
}

//------------------------------------------------------------------------------
// Initialization.
//------------------------------------------------------------------------------
FILE* GeomPolygon::mDevnull = NULL;

}


