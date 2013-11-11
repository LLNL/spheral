//---------------------------------Spheral++----------------------------------//
// GeomPolyhedron -- Geometric polyhedron class.
//
// A 2-D structure representing a polyhedron as a collection of GeomFacets.
//
// Created by JMO, Fri Jan 29 14:36:01 PST 2010
//----------------------------------------------------------------------------//
#include <algorithm>
#include <numeric>
#include <map>
#include <limits>

#include "boost/foreach.hpp"

#include "GeomPolyhedron.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/comparisons.hh"
#include "Utilities/PairComparisons.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/pointInPolyhedron.hh"

extern "C" {
#include "libqhull/qhull_a.h"
}

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
GeomPolyhedron::
GeomPolyhedron():
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
GeomPolyhedron::
GeomPolyhedron(const vector<GeomPolyhedron::Vector>& points):
  mVertices(),
  mFacets(),
  mConvex(true) {
  if (mDevnull == NULL) mDevnull = fopen("/dev/null", "w");

  if (points.size() > 0) {
    REQUIRE(points.size() > 3);

    // Find the appropriate renormalization so that we can give Qhull points
    // in a unit box.  Qhull just seems to work better this way.
    Vector xmin, xmax;
    boundingBox(points, xmin, xmax);
    const double fscale = (xmax - xmin).maxElement();
    CHECK(fscale > 0.0);

    // Copy the point coordinates to a Qhull point array.
    std::vector<coordT> points_qhull;
    points_qhull.reserve(3 * points.size());
    BOOST_FOREACH(Vector vec, points) {
      points_qhull.push_back((vec.x() - xmin.x())/fscale);
      points_qhull.push_back((vec.y() - xmin.y())/fscale);
      points_qhull.push_back((vec.z() - xmin.z())/fscale);
    }
    CHECK(points_qhull.size() == 3*points.size());

    // Call Qhull to generate the hull (C interface).
    boolT ismalloc = False;   /* True if qhull should free points in qh_freeqhull() or reallocation */
    char flags[250];          /* option flags for qhull, see qh_opt.htm */
//     FILE *outfile= NULL;      /* output from qh_produce_output()
//                                  use NULL to skip qh_produce_output() */
//     FILE *errfile= stderr;    /* error messages from qhull code */
//     FILE *errfile= fopen("/dev/null", "w");    /* error messages from qhull code */
    int exitcode;             /* 0 if no error from qhull */
    facetT *facet;            /* set by FORALLfacets */
    int curlong, totlong;     /* memory remaining after qh_memfreeshort */
    sprintf(flags, "qhull s"); // Tcv");
    const int exitcode_qhull = qh_new_qhull(3, points.size(), &points_qhull.front(), ismalloc, flags, mDevnull, mDevnull);

//     if (exitcode_qhull != 0) {
//       // Something didn't work, so spit out the points for diagnostics.
//       vector<pair<double, Vector> > sorted_points;
//       for (vector<Vector>::const_iterator itr = points.begin();
//            itr != points.end();
//            ++itr) sorted_points.push_back(make_pair(itr->magnitude(), *itr));
//       sort(sorted_points.begin(), sorted_points.end(), ComparePairByFirstElement<pair<double, Vector> >());
//       for (size_t k = 0; k != sorted_points.size(); ++k) cerr << "    -----> " << sorted_points[k].first << " " << sorted_points[k].second << endl;

//       // Emit the error message by calling qhull again.
//       FILE *errfile= stderr;    /* error messages from qhull code */
//       const int exitcode_qhull = qh_new_qhull(3, points.size(), &points_qhull.front(), ismalloc, flags, errfile, errfile);
//     }

    // VERIFY2(exitcode_qhull == 0,
    //         "Qhull emitted an error code while generating GeomPolyhedron");

    // If Qhull failed, we fall back on constructing the bounding box.
    if (exitcode_qhull != 0) {
      mVertices.push_back(Vector(xmax.x(), xmin.y(), xmin.z()));
      mVertices.push_back(Vector(xmax.x(), xmax.y(), xmin.z()));
      mVertices.push_back(Vector(xmin.x(), xmax.y(), xmin.z()));
      mVertices.push_back(Vector(xmin.x(), xmin.y(), xmin.z()));
      mVertices.push_back(Vector(xmax.x(), xmin.y(), xmax.z()));
      mVertices.push_back(Vector(xmax.x(), xmax.y(), xmax.z()));
      mVertices.push_back(Vector(xmin.x(), xmax.y(), xmax.z()));
      mVertices.push_back(Vector(xmin.x(), xmin.y(), xmax.z()));
      vector<unsigned> ipoints0, ipoints1, ipoints2, ipoints3, ipoints4, ipoints5;
      ipoints0.push_back(0); ipoints0.push_back(3); ipoints0.push_back(2); ipoints0.push_back(1);
      ipoints1.push_back(0); ipoints1.push_back(1); ipoints1.push_back(5); ipoints1.push_back(4);
      ipoints2.push_back(1); ipoints2.push_back(2); ipoints2.push_back(6); ipoints2.push_back(5);
      ipoints3.push_back(2); ipoints3.push_back(3); ipoints3.push_back(7); ipoints3.push_back(6);
      ipoints4.push_back(0); ipoints4.push_back(4); ipoints4.push_back(7); ipoints4.push_back(3);
      ipoints5.push_back(4); ipoints5.push_back(5); ipoints5.push_back(6); ipoints5.push_back(7);
      mFacets.push_back(Facet(mVertices, ipoints0, Vector(0, 0, -1)));
      mFacets.push_back(Facet(mVertices, ipoints1, Vector(1, 0, 0)));
      mFacets.push_back(Facet(mVertices, ipoints2, Vector(0, 1, 0)));
      mFacets.push_back(Facet(mVertices, ipoints3, Vector(-1, 0, 0)));
      mFacets.push_back(Facet(mVertices, ipoints4, Vector(0, -1, 0)));
      mFacets.push_back(Facet(mVertices, ipoints5, Vector(0, 0, 1)));

    } else {

      // Copy Qhull's vertex information.
      vertexT *vertex, **vertexp;
      map<unsigned, unsigned> vertexIDmap;
      {
        unsigned i = 0;
        FORALLvertices {
          mVertices.push_back(fscale*Vector(vertex->point[0], vertex->point[1], vertex->point[2]) + xmin);
          vertexIDmap[vertex->id] = i;
          ++i;
        }
      }
      CHECK(vertexIDmap.size() == mVertices.size());

      // Copy Qhull's facet information to our own format.
      vector<unsigned> vertex_indices;
      FORALLfacets {
        vertex_indices = vector<unsigned>();
        Vector fc;
        FOREACHvertex_(facet->vertices) {
          CHECK(vertexIDmap.find(vertex->id) != vertexIDmap.end());
          const unsigned vertex_i = vertexIDmap[vertex->id];
          CHECK(vertex_i < mVertices.size());
          vertex_indices.push_back(vertex_i);
          fc += mVertices[vertex_indices.back()];
        }
        CHECK(vertex_indices.size() >= 3);
        fc /= vertex_indices.size();
        const Vector normal(facet->normal[0], facet->normal[1], facet->normal[2]);

        // Enforce counter-clockwise ordering (viewed from outside) for the facet vertices.
        CHECK(vertex_indices.size() > 2);
        sort(vertex_indices.begin() + 1, vertex_indices.end(), CounterClockwiseComparator<Vector, vector<Vector> >(mVertices, mVertices[vertex_indices[0]], normal));

//       // If this is a triangular facet, we enforce counter-clockwise ordering (viewed from the outside) for the facet vertices.
//       if (vertex_indices.size() == 3) {
//         if ((mVertices[vertex_indices[1]] - fc).cross(mVertices[vertex_indices[0]] - fc).dot(normal) < 0.0) swap(vertex_indices[0], vertex_indices[1]);
//         if ((mVertices[vertex_indices[2]] - fc).cross(mVertices[vertex_indices[0]] - fc).dot(normal) < 0.0) swap(vertex_indices[0], vertex_indices[2]);
//         if ((mVertices[vertex_indices[2]] - fc).cross(mVertices[vertex_indices[1]] - fc).dot(normal) < 0.0) swap(vertex_indices[1], vertex_indices[2]);
//         CHECK((mVertices[vertex_indices[2]] - mVertices[vertex_indices[0]]).cross(mVertices[vertex_indices[1]] - mVertices[vertex_indices[0]]).dot(normal) > 0.0);
//       }

        mFacets.push_back(Facet(mVertices,
                                vertex_indices,
                                normal));
      }
    }

    // Free Qhull's resources.
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // Fill in our bounding box.
    setBoundingBox();

    // Post-conditions.
    BEGIN_CONTRACT_SCOPE;
    {
      // All normals should be outward facing.
      Vector centroid;
      for (vector<Vector>::const_iterator itr = mVertices.begin();
           itr != mVertices.end();
           ++itr) centroid += *itr;
      centroid /= mVertices.size();
      BOOST_FOREACH(const Facet& facet, mFacets) ENSURE2((facet.position() - centroid).dot(facet.normal()) >= 0.0,
                                                         "Inward normal? " << (facet.position() - centroid).dot(facet.normal()) << " " << facet.position() << " " << centroid << " " << facet.normal());

      // We had better be convex if built from a convex hull.
      ENSURE(convex());

//       // Ensure the seed points are contained.
//       for (vector<Vector>::const_iterator itr = points.begin();
//            itr != points.end();
//            ++itr) ENSURE(convexContains(*itr));
    }
    END_CONTRACT_SCOPE;

  }
}

//------------------------------------------------------------------------------
// Construct given the positions and facet indices.
// We assume here that the nodes for each facet are arranged correctly to
// create outward pointing normals.
//------------------------------------------------------------------------------
GeomPolyhedron::
GeomPolyhedron(const vector<GeomPolyhedron::Vector>& points,
               const vector<vector<unsigned> >& facetIndices):
  mVertices(points),
  mFacets(),
  mConvex(false) {

  unsigned i;
  vector<unsigned> indices;
  Vector a, b, normal, centroid;
  mFacets.reserve(facetIndices.size());

  // Construct the facets.
  BOOST_FOREACH(const vector<unsigned>& indices, facetIndices) {
    VERIFY2(indices.size() > 2, "Need at least two points per facet.");
    VERIFY2(*max_element(indices.begin(), indices.end()) < points.size(),
            "Bad vertex index for facet.");

    // Pick three points and construct an arbitrary normal.
    centroid.Zero();
    BOOST_FOREACH(i, indices) centroid += mVertices[i];
    centroid /= indices.size();
    a = mVertices[indices[0]] - centroid;
    b = mVertices[indices[1]] - centroid;
    normal = a.cross(b);
    VERIFY(normal.magnitude2() > 0.0);
    normal = normal.unitVector();

    mFacets.push_back(Facet(mVertices, indices, normal));
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
GeomPolyhedron::
GeomPolyhedron(const GeomPolyhedron& rhs):
  mVertices(rhs.mVertices),
  mFacets(),
  mXmin(rhs.mXmin),
  mXmax(rhs.mXmax),
  mConvex(rhs.mConvex) {
  mFacets.reserve(rhs.mFacets.size());
  BOOST_FOREACH(const Facet& facet, rhs.mFacets) mFacets.push_back(Facet(mVertices,
                                                                         facet.ipoints(),
                                                                         facet.normal()));
  ENSURE(mFacets.size() == rhs.mFacets.size());
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
GeomPolyhedron&
GeomPolyhedron::
operator=(const GeomPolyhedron& rhs) {
  if (this != &rhs) {
    mVertices = rhs.mVertices;
    mFacets = vector<Facet>();
    mFacets.reserve(rhs.mFacets.size());
    BOOST_FOREACH(const Facet& facet, rhs.mFacets) mFacets.push_back(Facet(mVertices,
                                                                           facet.ipoints(),
                                                                           facet.normal()));
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
GeomPolyhedron::
~GeomPolyhedron() {
}

//------------------------------------------------------------------------------
// Test if the given point is on the interior of the polyhedron.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
contains(const GeomPolyhedron::Vector& point,
         const bool countBoundary,
         const double tol) const {
  if (mConvex) {
    return this->convexContains(point, countBoundary, tol);
  } else {
    return pointInPolyhedron(point, *this, countBoundary, tol);
  }
}

//------------------------------------------------------------------------------
// Test if the given point is on the interior of the polyhedron.
// This method only works for convex polyhedra!
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
convexContains(const GeomPolyhedron::Vector& point,
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
// Test if we intersect the given polyhedron.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
intersect(const GeomPolyhedron& rhs) const {
  if (not testBoxIntersection(mXmin, mXmax, rhs.mXmin, rhs.mXmax)) return false;
  Vector vec;
  BOOST_FOREACH(vec, mVertices) {
    if (rhs.contains(vec)) return true;
  }
  BOOST_FOREACH(vec, rhs.mVertices) {
    if (this->contains(vec)) return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// Test if we intersect the given polyhedron assuming both are convex.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
convexIntersect(const GeomPolyhedron& rhs) const {
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
// Compute the centroid.
//------------------------------------------------------------------------------
GeomPolyhedron::Vector
GeomPolyhedron::
centroid() const {
  const unsigned n = mVertices.size();
  Vector result;
  if (n == 0) return result;
  CHECK(n >= 4);
  const Vector c0 = std::accumulate(mVertices.begin(), mVertices.end(), Vector())/n;

  // Specialize for tets, which are easy!
  if (n == 4) return c0;

  // Walk the facets.  Note we assume all facets are planar here!
  unsigned i, j;
  double vol, volsum = 0.0;
  Vector fc;
  BOOST_FOREACH(const Facet& facet, mFacets) {
    fc = facet.position();
    vol = facet.area() * facet.normal().dot(facet.point(0) - c0);  // Should be one third of this, but will cancel.
    volsum += vol;
    result += vol * (0.25*c0 + 0.75*fc);
  }
  CHECK(volsum > 0.0);
  result /= volsum;
  return result;
}

//------------------------------------------------------------------------------
// Return the edges of the polyhedron as a set of integer pairs for the vertices.
//------------------------------------------------------------------------------
vector<pair<unsigned, unsigned> >
GeomPolyhedron::
edges() const {
  vector<pair<unsigned, unsigned> > result;
  unsigned i, j, k, npoints;
  BOOST_FOREACH(const Facet& facet, mFacets) {
    const vector<unsigned>& ipoints = facet.ipoints();
    npoints = ipoints.size();
    for (k = 0; k != npoints; ++k) {
      i = ipoints[k];
      j = ipoints[(k + 1) % npoints];
      result.push_back(i < j ? make_pair(i, j) : make_pair(j, i));
    }
  }
  sort(result.begin(), result.end(), ComparePairs<pair<unsigned, unsigned> >());
  result.erase(unique(result.begin(), result.end()), result.end());
  return result;
}

//------------------------------------------------------------------------------
// Spit out an encoding of the facets as ordered vertex indices.
//------------------------------------------------------------------------------
vector<vector<unsigned> >
GeomPolyhedron::
facetVertices() const {
  vector<vector<unsigned> > result;
  if (mVertices.size() > 0) {
    BOOST_FOREACH(const Facet& facet, mFacets) {
      vector<unsigned> pts;
      copy(facet.ipoints().begin(), facet.ipoints().end(), back_inserter(pts));
      CHECK(pts.size() == facet.ipoints().size());
      result.push_back(pts);
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Spit out the facet normals.
//------------------------------------------------------------------------------
vector<GeomPolyhedron::Vector>
GeomPolyhedron::
facetNormals() const {
  vector<Vector> result;
  result.reserve(mFacets.size());
  BOOST_FOREACH(const Facet& facet, mFacets) result.push_back(facet.normal());
  ENSURE(result.size() == mFacets.size());
  return result;
}

//------------------------------------------------------------------------------
// Reconstruct the internal state given the set of vertices and the enocded 
// facets.
//------------------------------------------------------------------------------
void
GeomPolyhedron::
reconstruct(const vector<GeomPolyhedron::Vector>& vertices,
            const vector<vector<unsigned> >& facetVertices,
            const vector<GeomPolyhedron::Vector>& facetNormals) {
  const unsigned numFacets = facetVertices.size();
  REQUIRE(facetNormals.size() == numFacets);
  mVertices = vertices;
  mFacets = vector<Facet>();
  mFacets.reserve(numFacets);
  for (unsigned ifacet = 0; ifacet != numFacets; ++ifacet) {
    mFacets.push_back(Facet(mVertices,
                            facetVertices[ifacet],
                            facetNormals[ifacet]));
  }
  setBoundingBox();
  mConvex = this->convex();
  ENSURE(mFacets.size() == numFacets);
}

//------------------------------------------------------------------------------
// Compute the volume.
//------------------------------------------------------------------------------
double
GeomPolyhedron::
volume() const {
  double result = 0.0;
  const Vector c = centroid();
  BOOST_FOREACH(const Facet& facet, mFacets) {
    result += facet.area() * abs(facet.normal().dot(facet.point(0) - c));
  }
  ENSURE(result >= 0.0);
  return result/3.0;
}

//------------------------------------------------------------------------------
// Find the minimum distance to a point.
//------------------------------------------------------------------------------
double
GeomPolyhedron::
distance(const GeomPolyhedron::Vector& p) const {
  return (p - this->closestPoint(p)).magnitude();
}

//------------------------------------------------------------------------------
// Find the point in the polyhedron closest to the given point.
//------------------------------------------------------------------------------
GeomPolyhedron::Vector
GeomPolyhedron::
closestPoint(const GeomPolyhedron::Vector& p) const {
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
GeomPolyhedron::
operator==(const GeomPolyhedron& rhs) const {
  bool result = (mVertices == rhs.mVertices and
                 mFacets.size() == rhs.mFacets.size());
  int i = 0;
  while (result and i != mFacets.size()) {
    result = (mFacets[i] == rhs.mFacets[i]);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
operator!=(const GeomPolyhedron& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// Set the minimum and maximum extents of the polyhedron.
//------------------------------------------------------------------------------
void
GeomPolyhedron::
setBoundingBox() {
  boundingBox(mVertices, mXmin, mXmax);
}

//------------------------------------------------------------------------------
// Test if the polyhedron is convex.
//------------------------------------------------------------------------------
bool
GeomPolyhedron::
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
// Initialization.
//------------------------------------------------------------------------------
FILE* GeomPolyhedron::mDevnull = NULL;

}


