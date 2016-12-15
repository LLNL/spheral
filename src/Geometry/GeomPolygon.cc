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
#include <stdio.h>
#include <stdint.h>

#include "boost/foreach.hpp"

// #include "polytope/polytope.hh"
// #include "polytope/convexHull_2d.hh"

#include "GeomPolygon.hh"
#include "FacetedVolumeUtilities.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/pointInPolygon.hh"

//------------------------------------------------------------------------------
// It seems there is a missing specialization for abs(long unsigned int), so 
// fill it in.
// This is necessary for the collinear method below to compile.  It seems evil
// to insert something into namespace std:: like this, by the way.
//------------------------------------------------------------------------------
namespace std {
  // inline long unsigned int abs(long unsigned int x) { return x; }
  inline uint64_t          abs(uint64_t x)          { return x; }
}

namespace Spheral {

using namespace std;

//********************************************************************************
// The following anonymous stuff is lifted from the convex hull method I 
// implemented in polytope.
namespace {

struct KeyTraits {
  // typedef uint64_t Key;
  typedef int64_t Key;
  static const uint32_t numbits;
  static const uint32_t numbits1d;
  static const Key zero;
  static const Key one;
  static const Key two;
  static const Key maxKey1d;
  static const Key maxKey;
};

namespace geometry {

//------------------------------------------------------------------------------
// polytope 2D dot
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
RealType
dot(const RealType* a, const RealType* b) {
  return a[0]*b[0] + a[1]*b[1];
}

//------------------------------------------------------------------------------
// polytope 2D distance
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
RealType
distance(const RealType* a, const RealType* b) {
  const RealType dx = a[0] - b[0];
  const RealType dy = a[1] - b[1];
  return sqrt(dx*dx + dy*dy);
}

//------------------------------------------------------------------------------
// Determine if the given points are collinear to some accuracy.
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
collinear(const RealType* a, const RealType* b, const RealType* c, const RealType tol) {
  RealType ab[Dimension], ac[Dimension], abmag = 0.0, acmag = 0.0;
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] = b[j] - a[j];
    ac[j] = c[j] - a[j];
    abmag += ab[j]*ab[j];
    acmag += ac[j]*ac[j];
  }
  if (abmag < tol or acmag < tol) return true;
  abmag = sqrt(abmag);
  acmag = sqrt(acmag);
  for (unsigned j = 0; j != Dimension; ++j) {
    ab[j] /= abmag;
    ac[j] /= acmag;
  }
  return std::abs(std::abs(dot<Dimension, RealType>(ab, ac)) - 1.0) < tol;
}

}

//------------------------------------------------------------------------------
// A integer version of the simple 2D point.
//------------------------------------------------------------------------------
template<typename CoordType>
struct Point2 {
  CoordType x, y;
  unsigned index;
  Point2(): x(0), y(0), index(0) {}
  Point2(const CoordType& xi, const CoordType& yi, const unsigned i = 0): x(xi), y(yi), index(i) {}
  Point2& operator=(const Point2& rhs) { x = rhs.x; y = rhs.y; index = rhs.index; return *this; }
  bool operator==(const Point2& rhs) const { return (x == rhs.x and y == rhs.y); }
  bool operator!=(const Point2& rhs) const { return !(*this == rhs); }
  bool operator<(const Point2& rhs) const {
    return (x < rhs.x                ? true :
            x == rhs.x and y < rhs.y ? true :
            false);
  }
  template<typename RealType>
  Point2(const RealType& xi, const RealType& yi, const RealType& dx, const unsigned i = 0): 
    x(static_cast<CoordType>(xi/dx + 0.5)),
    y(static_cast<CoordType>(yi/dx + 0.5)),
    index(i) {}
  template<typename RealType>
  Point2(const RealType& xi, const RealType& yi, 
         const RealType& xlow, const RealType& ylow,
         const RealType& dx,
         const unsigned i = 0): 
    x(static_cast<CoordType>((xi - xlow)/dx + 0.5)),
    y(static_cast<CoordType>((yi - ylow)/dx + 0.5)),
    index(i) {}
  template<typename RealType> RealType realx(const RealType& xmin, const RealType& dx) const { return static_cast<RealType>(x*dx) + xmin; }
  template<typename RealType> RealType realy(const RealType& ymin, const RealType& dy) const { return static_cast<RealType>(y*dy) + ymin; }
  Point2& operator+=(const Point2& rhs) { x += rhs.x; y += rhs.y; return *this; }
  Point2& operator-=(const Point2& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
  Point2& operator*=(const CoordType& rhs) { x *= rhs; y *= rhs; return *this; }
  Point2& operator/=(const CoordType& rhs) { x /= rhs; y /= rhs; return *this; }
  Point2 operator+(const Point2& rhs) const { Point2 result(*this); result += rhs; return result; }
  Point2 operator-(const Point2& rhs) const { Point2 result(*this); result -= rhs; return result; }
  Point2 operator*(const CoordType& rhs) const { Point2 result(*this); result *= rhs; return result; }
  Point2 operator/(const CoordType& rhs) const { Point2 result(*this); result /= rhs; return result; }
  Point2 operator-() const { return Point2(-x, -y); }
  CoordType  operator[](const size_t i) const { CHECK(i < 2); return *(&x + i); }
  CoordType& operator[](const size_t i)       { CHECK(i < 2); return *(&x + i); }
};

//------------------------------------------------------------------------------
// A fuzzy comparison operator for our quantized Point2 type.
//------------------------------------------------------------------------------
template<typename UintType>
struct FuzzyPoint2LessThan {
  UintType fuzz;
  FuzzyPoint2LessThan(const UintType ifuzz = 1): fuzz(ifuzz) {}
  bool operator()(const Point2<UintType>& p1, const Point2<UintType>& p2) {
    return (p1.x + fuzz < p2.x        ? true :
            p1.x        > p2.x + fuzz ? false :
            p1.y + fuzz < p2.y        ? true :
            p1.y        > p2.y + fuzz ? false :
            false);
  }
  bool operator()(const std::pair<Point2<UintType>, unsigned>& p1,
                  const std::pair<Point2<UintType>, unsigned>& p2) {
    return operator()(p1.first, p2.first);
  }
};

//------------------------------------------------------------------------------
// sign of the Z coordinate of cross product : (p2 - p1)x(p3 - p1).
//------------------------------------------------------------------------------
template<typename RealType>
int zcross_sign(const Point2<RealType>& p1, const Point2<RealType>& p2, const Point2<RealType>& p3) {
//   double scale = 1.0/max(RealType(1), max(p1.x, max(p1.y, max(p2.x, max(p2.y, max(p3.x, p3.y))))));
  const double ztest = 
    (double(p2.x) - double(p1.x))*(double(p3.y) - double(p1.y)) -
    (double(p2.y) - double(p1.y))*(double(p3.x) - double(p1.x));
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
                         0);
  // return (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
}

//------------------------------------------------------------------------------
// Comparator to compare std::pair's by their first element.
//------------------------------------------------------------------------------
template<typename T1, typename T2>
struct ComparePairByFirstElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.first < rhs.first;
  }
};

//------------------------------------------------------------------------------
// The method itself.
//
// NOTE: The convex hull can be of dimension smaller than 2D. Lower
//       dimensionality is stored in the structure of the PLC facets
//             1D - Collinear points - A single length-2 facet with indices 
//                                     pointing to the smallest and largest
//                                     points in the sorted point hash
//             0D - Single point     - A single length-2 facet with the same
//                                     index (0) in both positions 
//------------------------------------------------------------------------------
template<typename RealType>
std::vector<std::vector<int> >
convexHull_2d(const std::vector<RealType>& points,
              const RealType* low,
              const RealType& dx) {
  typedef KeyTraits::Key CoordHash;
  typedef Point2<CoordHash> PointHash;
  // typedef polytope::DimensionTraits<2, RealType>::CoordHash CoordHash;
  // typedef polytope::DimensionTraits<2, RealType>::IntPoint PointHash;

  CHECK(!points.empty());
  CHECK(points.size() % 2 == 0);
  const unsigned n = points.size() / 2;
  vector<vector<int> > plc;
  int i, j, k, t;
  
  // If there's only one or two points, we're done: that's the whole hull
  if (n == 1 or n == 2) {
    plc.resize(1, std::vector<int>(2));
    plc[0][0] = 0;
    plc[0][1] = (n == 1) ? 0 : 1;
    return plc;
  }
  
  // Start by finding a point distinct from point 0.
  j = 1;
  while (j != n and geometry::distance<2, RealType>(&points[0], &points[2*j]) < dx) ++j;
  if (j == n - 1) {
    // There are only 2 distinct positions!
    plc.resize(1, std::vector<int>(2));
    plc[0][0] = 0;
    plc[0][1] = j;
    return plc;
  } else if (j == n) {
    // Good god, there are no distinct points!
    plc.resize(1, std::vector<int>(2));
    plc[0][0] = 0;
    plc[0][1] = 0;
    return plc;
  }

  // Check if the input points are collinear.
  bool collinear = true;
  CHECK(n > 2);
  i = 2;
  while (collinear and i != n) {
    collinear = geometry::collinear<2,RealType>(&points[0], &points[2*j], &points[2*i], dx);
    ++i;
  }
  
  // Hash the input points and sort them by x coordinate, remembering their original indices
  // in the input set.  We also ensure that only unique (using a fuzzy comparison) points
  // are inserted here, since duplicates mess up the hull calculation.
  const RealType& xmin = low[0];
  const RealType& ymin = low[1];
  std::set<std::pair<PointHash, unsigned>, FuzzyPoint2LessThan<CoordHash> > uniquePoints;
  for (i = 0; i != n; ++i) {
    uniquePoints.insert(std::make_pair(PointHash(CoordHash((points[2*i]     - xmin)/dx + 0.5),
                                                 CoordHash((points[2*i + 1] - ymin)/dx + 0.5)),
                                       i));
  }
  std::vector<std::pair<PointHash, unsigned> > sortedPoints(uniquePoints.begin(), uniquePoints.end());
  std::sort(sortedPoints.begin(), sortedPoints.end());

  // If the points are collinear, we can save a lot of work
  if (collinear) {
    plc.resize(1, std::vector<int>(2));
    plc[0][0] = sortedPoints.front().second;
    plc[0][1] = sortedPoints.back().second;
  }
  else {
    // Prepare the result.
    const unsigned nunique = sortedPoints.size();
    std::vector<int> result(2*nunique);
    
    // Build the lower hull.
    for (i = 0, k = 0; i < nunique; i++) {
      while (k >= 2 and
             zcross_sign(sortedPoints[result[k - 2]].first, sortedPoints[result[k - 1]].first, sortedPoints[i].first) <= 0) k--;
      result[k++] = i;
    }
    
    // Build the upper hull.
    for (i = nunique - 2, t = k + 1; i >= 0; i--) {
      while (k >= t and
             zcross_sign(sortedPoints[result[k - 2]].first, sortedPoints[result[k - 1]].first, sortedPoints[i].first) <= 0) k--;
      result[k++] = i;
    }
    // if (!(k >= 4)) {
    //   std::cerr << "Blago!  " << n << " " << nunique << " " << k << std::endl;
    //   std::cerr << "Unique:" << std::endl;
    //   for (unsigned i = 0; i != nunique; ++i) std::cerr << "  --> " << sortedPoints[i].first << std::endl;
    //   std::cerr << "Input:" << std::endl;
    //   for (unsigned i = 0; i != n; ++i) std::cerr << "  --> " << points[2*i] << " " << points[2*i+1] << std::endl;
    // }
    CHECK(k >= 4);
    CHECK(result.front() == result.back());
    
    // Translate our sorted information to a PLC based on the input point ordering and we're done.
    for (i = 0; i != k - 1; ++i) {
      j = (i + 1) % k;
      plc.push_back(std::vector<int>());
      plc.back().push_back(sortedPoints[result[i]].second);
      plc.back().push_back(sortedPoints[result[j]].second);
    }
    CHECK(plc.size() == k - 1);
  }
  return plc;
}

} // end anonymous namespace
//********************************************************************************

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
GeomPolygon::
GeomPolygon():
  mVertices(),
  mFacets(),
  mVertexUnitNorms(),
  mVertexFacetConnectivity(),
  mFacetFacetConnectivity(),
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
  mVertexUnitNorms(),
  mVertexFacetConnectivity(),
  mFacetFacetConnectivity(),
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
    double low[2] = {0.0, 0.0};
    // polytope::PLC<2, double> plc = polytope::convexHull_2d(points_polytope, &(*low.begin()), 1.0e-15);
    vector<vector<int> > plc = convexHull_2d(points_polytope, low, 1.0e-8);
    const unsigned numVertices = plc.size();
    CHECK(numVertices >= 3);

    // Extract the hull information back to our local convention.  We use the fact that
    // polytope's convex hull method sorts the vertices in counter-clockwise here.
    // Start with the vertices.
    mVertices.reserve(numVertices);
    int i, j;
    for (j = 0; j != numVertices; ++j) {
      CHECK(plc[j].size() == 2);
      i = plc[j][0];
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

    // Compute the ancillary geometry.
    GeometryUtilities::computeAncillaryGeometry(*this, mVertexFacetConnectivity, mFacetFacetConnectivity, mVertexUnitNorms);

    // Post-conditions.
    BEGIN_CONTRACT_SCOPE
    {
      // Ensure the facet node ordering is correct.
      CounterClockwiseComparator<Vector, vector<Vector> > nodeComparator(mVertices, this->centroid());
      BOOST_FOREACH(const Facet& facet, mFacets) ENSURE2(nodeComparator(facet.point1(), facet.point2()), *this);

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
        ENSURE(nodeComparator(i, j));
      }

      // We had better be convex if built from a convex hull.
      ENSURE(this->convex());

      // Ensure the seed points are contained.
      // Suspending this check for now as floating point accuracy occasionally misfires
      // this check.
      //      BOOST_FOREACH(vec, points) ENSURE(this->convexContains(vec));
    }
    END_CONTRACT_SCOPE
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

  // Compute the ancillary geometry.
  GeometryUtilities::computeAncillaryGeometry(*this, mVertexFacetConnectivity, mFacetFacetConnectivity, mVertexUnitNorms);
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
GeomPolygon::
GeomPolygon(const GeomPolygon& rhs):
  mVertices(rhs.mVertices),
  mFacets(),
  mVertexFacetConnectivity(rhs.mVertexFacetConnectivity),
  mFacetFacetConnectivity(rhs.mFacetFacetConnectivity),
  mVertexUnitNorms(rhs.mVertexUnitNorms),
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
    BOOST_FOREACH(const Facet& facet, rhs.mFacets) mFacets.push_back(Facet(mVertices,
                                                                           facet.ipoint1(),
                                                                           facet.ipoint2()));
    mVertexFacetConnectivity = rhs.mVertexFacetConnectivity;
    mFacetFacetConnectivity = rhs.mVertexFacetConnectivity;
    mVertexUnitNorms = rhs.mVertexUnitNorms;
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
    return this->convexContains(point, countBoundary, tol);
  } else {
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
      ++facetItr;
    }
  } else {
    while (facetItr != mFacets.end() and result) {
      result = (facetItr->compare(point, tol) < 0);
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
// Test if we intersect the given box.
//------------------------------------------------------------------------------
bool
GeomPolygon::
intersect(const std::pair<Vector, Vector>& rhs) const {
  if (not testBoxIntersection(mXmin, mXmax, rhs.first, rhs.second)) return false;
  
  // Build a GeompPolygon representation of the box and use our generic intersection
  // method.
  vector<Vector> verts(4);
  verts[0] = rhs.first; 
  verts[1] = Vector(rhs.second.x(), rhs.first.y());
  verts[2] = rhs.second;
  verts[3] = Vector(rhs.first.x(), rhs.second.y());
  vector<vector<unsigned> > facets(4);
  facets[0].push_back(0); facets[0].push_back(1);
  facets[1].push_back(1); facets[1].push_back(2);
  facets[2].push_back(2); facets[2].push_back(3);
  facets[3].push_back(3); facets[3].push_back(0);
  GeomPolygon other(verts, facets);
  return this->intersect(other);
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
  GeometryUtilities::computeAncillaryGeometry(*this, mVertexFacetConnectivity, mFacetFacetConnectivity, mVertexUnitNorms);
  ENSURE(mFacets.size() == facetVertices.size());
  ENSURE(mFacetFacetConnectivity.size() == mFacets.size());
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
// Find the facet closest to the given point.
//------------------------------------------------------------------------------
unsigned
GeomPolygon::
closestFacet(const GeomPolygon::Vector& p) const {
  unsigned result = 0;
  double r2, minr2 = numeric_limits<double>::max();
  Vector thpt;
  for (unsigned i = 0; i != mFacets.size(); ++i) {
    thpt = mFacets[i].closestPoint(p);
    r2 = (thpt - p).magnitude2();
    if (r2 < minr2) {
      result = i;
      minr2 = r2;
    }
  }
  ENSURE(result < mFacets.size());
  return result;
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
  const Facet& f = mFacets[this->closestFacet(p)];
  return f.closestPoint(p);
}

//------------------------------------------------------------------------------
// += Vector, shift polygon in space
//------------------------------------------------------------------------------
GeomPolygon&
GeomPolygon::
operator+=(const GeomPolygon::Vector& rhs) {
  for (auto& v: mVertices) v += rhs;
  mXmin += rhs;
  mXmax += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// -= Vector, shift polygon in space
//------------------------------------------------------------------------------
GeomPolygon&
GeomPolygon::
operator-=(const GeomPolygon::Vector& rhs) {
  (*this) += -rhs;
  return *this;
}

//------------------------------------------------------------------------------
// + Vector, return shifted polygon in space
//------------------------------------------------------------------------------
GeomPolygon
GeomPolygon::
operator+(const GeomPolygon::Vector& rhs) const {
  GeomPolygon result(*this);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// - Vector, return shifted polygon in space
//------------------------------------------------------------------------------
GeomPolygon
GeomPolygon::
operator-(const GeomPolygon::Vector& rhs) const {
  return (*this) + (-rhs);
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


