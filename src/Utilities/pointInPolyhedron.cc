//------------------------------------------------------------------------------
// pointInPolyhedron
//
// Test if a given point is in a polyhedron or not.
//------------------------------------------------------------------------------
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>

#include "pointInPolyhedron.hh"
#include "Geometry/Dimension.hh"
#include "testBoxIntersection.hh"
#include "pointOnPolyhedron.hh"
#include "segmentIntersectEdges.hh"
#include "pointInPolygon.hh"
#include "rotationMatrix.hh"
#include "lineSegmentIntersections.hh"

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

  // First reject any points outside the bounding box.
  if (not testPointInBox(p, polyhedron.xmin(), polyhedron.xmax(), tol)) return false;

  // Are we checking the surface?
  vector<Vector> vertices = polyhedron.vertices();
  if (pointOnPolyhedron(p, polyhedron, tol)) return countBoundary;

  // We're going to take O'Rourke's advice here and look for a ray that doesn't 
  // have any degeneracies -- i.e., it doesn't intersect any edges or vertices
  // of the polyhedron.
  const double length = (polyhedron.xmax() - polyhedron.xmin()).magnitude();
  Vector p1 = p + Vector(length, 0.0, 0.0);
  // {
  //   // Construct a random number generator.
  //   boost::mt19937 seed;
  //   boost::uniform_real<double> gm(-1.0, 1.0);
  //   boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > gen(seed, gm);

  //   const unsigned maxSearchIterations = 200;
  //   unsigned iter = 0;
  //   double theta, phi;
  //   while (++iter < maxSearchIterations and segmentIntersectEdges(p, p1, polyhedron))
  //     p1 = p + length * (Vector(gen(), gen(), gen()).unitVector());
  //   if (!(iter < maxSearchIterations)) {
  //     cerr << "Point:  Vector(" << p.x() << ", " << p.y() << ", " << p.z() << ")" << endl;
  //     cerr << "vertices:  [";
  //     for (unsigned i = 0; i != vertices.size(); ++i) cerr << "Vector(" << vertices[i].x() << ", " << vertices[i].y() << ", " << vertices[i].z() << "), ";
  //     cerr << "]" << endl;
  //   }
  //   CHECK2(iter < maxSearchIterations, "Failed to find a safe ray!  " << maxSearchIterations << " " << p1 - p << " " << p);
  // }
  // const Vector pp1 = p1 - p;
  //   cerr << "Chose ray:  " << pp1 << endl;

  // Walk the facets, and test each one.
  bool filter;
  unsigned numIntersections = 0;
  const vector<Facet>& facets = polyhedron.facets();
  bool facetTest;
  unsigned ifacet, npts, i, j;
  double px, py, pz, fxmin, fymin, fzmin, fxmax, fymax, fzmax;
  char code;
  Vector planeIntercept, normal;
  for (ifacet = 0; ifacet != facets.size(); ++ifacet) {
    const Facet& facet = facets[ifacet];
    const vector<unsigned>& ipts = facet.ipoints();
    npts = ipts.size();

    // Find the bounding box for the facet.
    fxmin =  1e50; fymin =  1e50; fzmin =  1e50;
    fxmax = -1e50; fymax = -1e50; fzmax = -1e50;
    for (i = 0; i != npts; ++i) {
      fxmin = min(fxmin, vertices[ipts[i]].x());
      fymin = min(fymin, vertices[ipts[i]].y());
      fzmin = min(fzmin, vertices[ipts[i]].z());
      fxmax = max(fxmax, vertices[ipts[i]].x());
      fymax = max(fymax, vertices[ipts[i]].y());
      fzmax = max(fzmax, vertices[ipts[i]].z());
    }

    // NOTE!  This culling stage only works so long as we are checking rays in the (1,0,0) direction
    // like currently hardwired.
    if (fxmax >= p.x() - tol and
        fymin <= p.y() + tol and fymax >= p.y() - tol and
        fzmin <= p.z() + tol and fzmax >= p.z() - tol) {

      normal = facet.normal().unitVector();
      //     cerr << "Point-facet distance: " << pointPlaneDistance(p, facet.position(), normal) << endl;

      // Find the intersection of the line corresponding to our ray and the facet plane.
      code = segmentPlaneIntersection(p, p1, facet.position(), normal, planeIntercept);
      //     cerr << "Code : " << code << " " << (planeIntercept - facet.position()).dot(normal) << " " << (planeIntercept - p).dot(pp1) << endl;
      if (code != '0') { // and (planeIntercept - p).dot(pp1) >= 0.0) {
        facetTest = false;
        px = planeIntercept.x();
        py = planeIntercept.y();
        pz = planeIntercept.z();

        // The plane intersection is in the positive direction of the ray, so we 
        // project the facet and intercept point to one of the primary planes and
        // do a 2-D polygon interior test.
        if (abs(normal.x()) > 0.5) {

          // x plane -- use (y,z) coordinates.
          if (py >= fymin and py <= fymax and pz >= fzmin and pz <= fzmax) {
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].z() > pz) != (vertices[ipts[j]].z() > pz)) &&
                   (py < (vertices[ipts[j]].y() - vertices[ipts[i]].y()) * (pz - vertices[ipts[i]].z()) / (vertices[ipts[j]].z() - vertices[ipts[i]].z()) + vertices[ipts[i]].y()) )
                facetTest = not facetTest;
            }
          }

          //         cerr << "x test: " << facetTest << endl;
          //         vector<Dim<2>::Vector> v2d;
          //         for (i = 0; i != npts; ++i) v2d.push_back(Dim<2>::Vector(vertices[ipts[i]].y(), vertices[ipts[i]].z()));
          //         cerr << "Polygonal test:  " << pointInPolygon(Dim<2>::Vector(py, pz), v2d, false) << " " << pointInPolygon(Dim<2>::Vector(py, pz), v2d, true) << endl;

        } else if (abs(normal.y()) > 0.5) {

          // y plane -- use (x,z) coordinates.
          if (px >= fxmin and px <= fxmax and pz >= fzmin and pz <= fzmax) {
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].z() > pz) != (vertices[ipts[j]].z() > pz)) &&
                   (px < (vertices[ipts[j]].x() - vertices[ipts[i]].x()) * (pz - vertices[ipts[i]].z()) / (vertices[ipts[j]].z() - vertices[ipts[i]].z()) + vertices[ipts[i]].x()) )
                facetTest = not facetTest;
            }
          }

          //         cerr << "y test: " << facetTest << endl;
          //         vector<Dim<2>::Vector> v2d;
          //         for (i = 0; i != npts; ++i) v2d.push_back(Dim<2>::Vector(vertices[ipts[i]].x(), vertices[ipts[i]].z()));
          //         cerr << "Polygonal test:  " << pointInPolygon(Dim<2>::Vector(px, pz), v2d, false) << " " << pointInPolygon(Dim<2>::Vector(px, pz), v2d, true) << endl;

        } else {
          CHECK(abs(normal.z()) > 0.5);

          // z plane -- use (x,y) coordinate.
          if (px >= fxmin and px <= fxmax and py >= fymin and py <= fymax) {
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].y() > py) != (vertices[ipts[j]].y() > py)) &&
                   (px < (vertices[ipts[j]].x() - vertices[ipts[i]].x()) * (py - vertices[ipts[i]].y()) / (vertices[ipts[j]].y() - vertices[ipts[i]].y()) + vertices[ipts[i]].x()) )
                facetTest = not facetTest;
            }
          }

          //         cerr << "z test: " << facetTest << endl;
          //         vector<Dim<2>::Vector> v2d;
          //         for (i = 0; i != npts; ++i) v2d.push_back(Dim<2>::Vector(vertices[ipts[i]].x(), vertices[ipts[i]].y()));
          //         cerr << "Polygonal test:  " << pointInPolygon(Dim<2>::Vector(px, py), v2d, false) << " " << pointInPolygon(Dim<2>::Vector(px, py), v2d, true) << endl;

        }
        if (facetTest) ++numIntersections;

        //       cerr << "Code: " << code << " " << p << " " << planeIntercept << " " << (planeIntercept - p) << " " << (planeIntercept - p).dot(pp1) << endl;
        //       if (facetTest) {
        //         cerr << "Facet:  " << numIntersections << " " << facetTest << " " << facet.position() << " " << facet.normal() << endl;
        //         cerr << "        " << pointPlaneDistance(p, facet.position(), facet.normal()) << " : " << endl;
        //       }
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
//   }

  return not (numIntersections % 2 == 0);
}

}
