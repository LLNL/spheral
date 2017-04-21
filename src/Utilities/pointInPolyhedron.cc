//------------------------------------------------------------------------------
// pointInPolyhedron
//
// Test if a given point is in a polyhedron or not.
//------------------------------------------------------------------------------
#include <random>
#include <algorithm>

#include "pointInPolyhedron.hh"
#include "Geometry/Dimension.hh"
#include "testBoxIntersection.hh"
#include "pointOnPolyhedron.hh"
#include "segmentIntersectEdges.hh"
#include "lineSegmentIntersections.hh"
#include "rotationMatrix.hh"
#include "boost/random/uniform_real_distribution.hpp"

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// The method in question -- test if a point is on the interior of a polyhedron.
//------------------------------------------------------------------------------
bool pointInPolyhedron(const Dim<3>::Vector& p,
                       const Dim<3>::FacetedVolume& polyhedron,
                       const bool countBoundary,
                       const double tol) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::FacetedVolume Polyhedron;
  typedef Polyhedron::Facet Facet;
  typedef Dim<2>::Vector Vector2d;

  const unsigned numTrials = 1;
  const vector<Vector>& vertices = polyhedron.vertices();

  // First reject any points outside the bounding box.
  if (not testPointInBox(p, polyhedron.xmin(), polyhedron.xmax(), tol)) return false;

  // Are we checking the surface?
  if (countBoundary and pointOnPolyhedron(p, polyhedron, tol)) return true;

  // Build our ray for projecting through the polyhedron.
  const double length = std::max((polyhedron.xmax() - polyhedron.xmin()).magnitude(),
                                 std::max((p - polyhedron.xmin()).magnitude(),
                                          (p - polyhedron.xmax()).magnitude())) * 10.0;
  // Vector p1 = Vector(length, 0.0, 0.0);
  // Vector p2 = Vector(0.0, length, 0.0);
  // Vector p3 = Vector(0.0, 0.0, length);

  // Pick a random rotation to apply to these rays.
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  boost::random::uniform_real_distribution<> dis(-1.0, 1.0);
  // const Tensor RT = rotationMatrix(Vector(dis(gen), dis(gen), dis(gen)).unitVector());
  // p1 = p + RT*p1;
  // p2 = p + RT*p2;
  // p3 = p + RT*p3;

  // We're going to take O'Rourke's advice here and look for a ray that doesn't 
  // have any degeneracies -- i.e., it doesn't intersect any edges or vertices
  // of the polyhedron.
  std::vector<Vector> pends;
  while (pends.size() < numTrials) {
    auto iter = 0;
    Tensor RT = rotationMatrix(Vector(dis(gen), dis(gen), dis(gen)).unitVector());
    Vector p1 = p + RT*Vector(length, 0.0, 0.0);
    while (iter < 200 and segmentIntersectEdges(p1, p1, polyhedron)) {
      RT = rotationMatrix(Vector(dis(gen), dis(gen), dis(gen)).unitVector());
      p1 = p + RT*Vector(length, 0.0, 0.0);
    }
    VERIFY2(iter < 200, "Failed to cast a non-degenerate ray through polyhedron.");
    pends.push_back(p1);
  }

  // Walk the facets, and test each one.
  vector<unsigned> numIntersections(numTrials, 0);
  const vector<Facet>& facets = polyhedron.facets();
  bool filter, facetTest, boundTest;
  unsigned ifacet, npts, i, j;
  double px, py, pz, fxmin, fymin, fzmin, fxmax, fymax, fzmax;
  char code;
  Vector planeIntercept, normal;
  
  for (ifacet = 0; ifacet != facets.size(); ++ifacet) {
    const Facet& facet = facets[ifacet];
    const vector<unsigned>& ipts = facet.ipoints();
    npts = ipts.size();

    for (auto itrial = 0; itrial < numTrials; ++itrial) {
      Vector ptest = pends[itrial];

      filter = false;
      i = 0;
      while (i != npts and not filter) {
        filter = ((vertices[ipts[i++]] - p).dot(ptest - p) >= 0.0);
      }

      if (filter) {

        normal = facet.normal().unitVector();

        // Find the intersection of the line corresponding to our ray and the facet plane.
        code = segmentPlaneIntersection(p, ptest, facet.position(), normal, planeIntercept, tol);
        if (code != '0') { // and (planeIntercept - p).dot(pp1) >= 0.0) {

          // Test the interior of the facet.
          facetTest = false;
          px = planeIntercept.x();
          py = planeIntercept.y();
          pz = planeIntercept.z();

          // The plane intersection is in the positive direction of the ray, so we 
          // project the facet and intercept point to one of the primary planes and
          // do a 2-D polygon interior test.
          if (abs(normal.x()) > 0.1) {

            // x plane -- use (y,z) coordinates.
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].z() > pz) != (vertices[ipts[j]].z() > pz)) &&
                   (py < (vertices[ipts[j]].y() - vertices[ipts[i]].y()) * (pz - vertices[ipts[i]].z()) / (vertices[ipts[j]].z() - vertices[ipts[i]].z()) + vertices[ipts[i]].y()) )
                facetTest = not facetTest;
            }

          } else if (abs(normal.y()) > 0.1) {
       
            // y plane -- use (z,x) coordinates.
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].x() > px) != (vertices[ipts[j]].x() > px)) &&
                   (pz < (vertices[ipts[j]].z() - vertices[ipts[i]].z()) * (px - vertices[ipts[i]].x()) / (vertices[ipts[j]].x() - vertices[ipts[i]].x()) + vertices[ipts[i]].z()) )
                facetTest = not facetTest;
            }

          } else {
            CHECK(abs(normal.z()) > 0.1);

            // z plane -- use (x,y) coordinate.
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].y() > py) != (vertices[ipts[j]].y() > py)) &&
                   (px < (vertices[ipts[j]].x() - vertices[ipts[i]].x()) * (py - vertices[ipts[i]].y()) / (vertices[ipts[j]].y() - vertices[ipts[i]].y()) + vertices[ipts[i]].x()) )
                facetTest = not facetTest;
            }

          }
          if (facetTest) numIntersections[itrial] += 1;
        }
      }
    }
  }

  // If the ray passed through an even number of intersections the point is external to the 
  // polyhedron.
//   cerr << " ------> numIntersections:  " << numIntersections << endl
//        << "  Facet tests:  ";
//   {
//     int j = 0;
//     for (i = 0; i != facets.size(); ++i) {
//       j += facets[i].compare(p, tol);
//       cerr << " " << j;
//     }
//     cerr << endl
//          << testPointInBox(p, polyhedron.xmin(), polyhedron.xmax(), tol) << endl;
// }

  // VERIFY2((numIntersections1 % 2)  == (numIntersections2 % 2) and
  //         (numIntersections1 % 2)  == (numIntersections3 % 2),
  //         numIntersections1 << " " << numIntersections2 << " " << numIntersections3);
  // return ((numIntersections1 % 2 != 0) and (numIntersections2 % 2 != 0) and (numIntersections3 % 2 != 0));

  // Really all three of these measures should give us the same answer.  However, due to precision and such they can vary.  We punt here a bit and
  // take the best two out three as our answer.
  unsigned numtrue = 0, numfalse = 0;
  for (const auto trial: numIntersections) {
    if (trial % 2 == 0) {
      numfalse += 1;
    } else {
      numtrue += 1;
    }
  }
  return numtrue > numfalse;


  // return ((numIntersections1 % 2 != 0 ? 1 : 0) +
  //         (numIntersections2 % 2 != 0 ? 1 : 0) +
  //         (numIntersections3 % 2 != 0 ? 1 : 0) >= 2);
}

}
