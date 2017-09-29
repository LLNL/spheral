//------------------------------------------------------------------------------
// pointInPolyhedron
//
// Test if a given point is in a polyhedron or not.
//------------------------------------------------------------------------------
#include <random>
#include <algorithm>

#include "pointInPolyhedron.hh"
#include "pointInPolygon.hh"
#include "testBoxIntersection.hh"
#include "pointOnPolyhedron.hh"
#include "segmentIntersectEdges.hh"
#include "lineSegmentIntersections.hh"
#include "rotationMatrix.hh"
#include "Geometry/Dimension.hh"

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

  const unsigned numTrials = 50;
  const auto& vertices = polyhedron.vertices();

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
  std::uniform_real_distribution<> dis(-1.0, 1.0);
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
    auto RT = rotationMatrix(Vector(dis(gen), dis(gen), dis(gen)).unitVector());
    auto p1 = p + RT*Vector(length, 0.0, 0.0);
    while (iter < 200 and segmentIntersectEdges(p, p1, polyhedron)) {
      RT = rotationMatrix(Vector(dis(gen), dis(gen), dis(gen)).unitVector());
      p1 = p + RT*Vector(length, 0.0, 0.0);
    }
    VERIFY2(iter < 200, "Failed to cast a non-degenerate ray through polyhedron.");
    pends.push_back(p1);
  }

  // Walk the facets, and test each one.
  vector<unsigned> numIntersections(numTrials, 0);
  const auto& facets = polyhedron.facets();
  bool xtest, ytest, ztest;
  unsigned i, j, npts;
  double px, py, pz;
  char code;
  Vector intercept1, intercept2, normal;
  
  for (const auto& facet: facets) {
    const auto& ipts = facet.ipoints();
    npts = ipts.size();

    for (auto itrial = 0; itrial < numTrials; ++itrial) {
      normal = facet.normal().unitVector();

      // Find the intersection (if any) of the line corresponding to our ray and the facet plane.
      code = segmentPlaneIntersection(p, pends[itrial], facet.position(), normal, intercept1, tol);
      CHECK(code != 'd');

      if (code == 'p') {

        // The line segement is in the plane of facet.
        if (pointInPolygon(p, vertices, ipts, normal, false, tol) or 
            pointInPolygon(pends[itrial], vertices, ipts, normal, false, tol)) {
          numIntersections[itrial] += 1;
        } else {
          i = 0;
          while (i < npts - 1 and 
                 segmentSegmentIntersection(p, pends[itrial], vertices[ipts[i]], vertices[ipts[i+1]], intercept1, intercept2, tol) == '0') ++i;
          if (i < npts - 1) numIntersections[itrial] += 1;
        }

      } else if (code == '1') { // and (intercept1 - p).dot(pp1) >= 0.0) {

        if (pointInPolygon(intercept1, vertices, ipts, normal, false, tol)) numIntersections[itrial] += 1;

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
    // if (Process::getRank() == 0) cout << "Trial : " << trial << " " << (trial % 2) << endl;
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
