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
  if (countBoundary and pointOnPolyhedron(p, polyhedron, tol)) return true;

  // Build our ray for projecting through the polyhedron.
  const double length = (polyhedron.xmax() - polyhedron.xmin()).magnitude();
  Vector p1 = p + Vector(length, 0.0, 0.0);
  Vector p2 = p + Vector(0.0, length, 0.0);
  Vector p3 = p + Vector(0.0, 0.0, length);

  // We're going to take O'Rourke's advice here and look for a ray that doesn't 
  // have any degeneracies -- i.e., it doesn't intersect any edges or vertices
  // of the polyhedron.
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
  unsigned numIntersections1 = 0, numIntersections2 = 0, numIntersections3 = 0;
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

    Vector ptest = p1;
    {
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
       
            // y plane -- use (x,z) coordinates.
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].z() > pz) != (vertices[ipts[j]].z() > pz)) &&
                   (px < (vertices[ipts[j]].x() - vertices[ipts[i]].x()) * (pz - vertices[ipts[i]].z()) / (vertices[ipts[j]].z() - vertices[ipts[i]].z()) + vertices[ipts[i]].x()) )
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
          if (facetTest) numIntersections1 += 1;
        }
      }
    }

    ptest = p2;
    {
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
       
            // y plane -- use (x,z) coordinates.
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].z() > pz) != (vertices[ipts[j]].z() > pz)) &&
                   (px < (vertices[ipts[j]].x() - vertices[ipts[i]].x()) * (pz - vertices[ipts[i]].z()) / (vertices[ipts[j]].z() - vertices[ipts[i]].z()) + vertices[ipts[i]].x()) )
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
          if (facetTest) numIntersections2 += 1;
        }
      }
    }

    ptest = p3;
    {
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
       
            // y plane -- use (x,z) coordinates.
            for (i = 0, j = npts - 1; i < npts; j = i++) {
              if ( ((vertices[ipts[i]].z() > pz) != (vertices[ipts[j]].z() > pz)) &&
                   (px < (vertices[ipts[j]].x() - vertices[ipts[i]].x()) * (pz - vertices[ipts[i]].z()) / (vertices[ipts[j]].z() - vertices[ipts[i]].z()) + vertices[ipts[i]].x()) )
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
          if (facetTest) numIntersections3 += 1;
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

  return ((numIntersections1 % 2 != 0) and (numIntersections2 % 2 != 0) and (numIntersections3 % 2 != 0));
}

}
