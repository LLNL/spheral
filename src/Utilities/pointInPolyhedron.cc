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

  const auto& vertices = polyhedron.vertices();
  const auto& facets = polyhedron.facets();
  const auto pz = p.z();
  const Vector2d pxy(p.x(), p.y()),
                 pzx(p.z(), p.x()),
                 pyz(p.y(), p.z());

  // First reject any points outside the bounding box.
  if (not testPointInBox(p, polyhedron.xmin(), polyhedron.xmax(), tol)) return false;

  // Are we checking the surface?
  if (countBoundary and polyhedron.distance(p) < tol) return true;
  // if (countBoundary and pointOnPolyhedron(p, polyhedron, tol)) return true;

  // Keep count of how mny facets we intersect.
  int numxy = 0, numzx = 0, numyz = 0;

  // Walk the facets.
  for (const auto& facet: facets) {
    const auto& ipts = facet.ipoints();
    const auto npts = ipts.size();

    // We only have to check if some portion of this facet is at greater z than the point in question.
    size_t i = 0;
    while (i < npts and vertices[ipts[i]].z() < pz) ++i;
    if (i < npts) {

      // Copy the vertices as XY pairs.
      vector<Vector2d> vxy, vzx, vyz;
      for (i = 0; i < npts; ++i) {
        vxy.push_back(Vector2d(vertices[ipts[i]].x(), vertices[ipts[i]].y()));
        vzx.push_back(Vector2d(vertices[ipts[i]].z(), vertices[ipts[i]].x()));
        vyz.push_back(Vector2d(vertices[ipts[i]].y(), vertices[ipts[i]].z()));
      }

      // Does a ray straight up the +z direction intersect this projected polygon?
      if (pointInPolygon(pxy, vxy, false, tol)) ++numxy;
      if (pointInPolygon(pzx, vzx, false, tol)) ++numzx;
      if (pointInPolygon(pyz, vyz, false, tol)) ++numyz;
    }
  }

  // If the ray passed through an even number of intersections the point is external to the 
  // polyhedron.
  const auto xytest = (numxy % 2) == 1,
             zxtest = (numzx % 2) == 1,
             yztest = (numyz % 2) == 1;
  if (xytest != zxtest or xytest != yztest or zxtest != yztest) {
    cerr << "Whoa : " << numxy << " " << numzx << " " << numyz << endl;
  }
  return zxtest;
  // return ((xytest and zxtest) or
  //         (zxtest and yztest) or
  //         (yztest and xytest));
}

}
