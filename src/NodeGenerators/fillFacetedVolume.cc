#include "fillFacetedVolume.hh"

namespace Spheral {

using namespace std;

//------------------------------------------------------------------------------
// Fill an outer bounding volume.
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  vector<Vector> result;
  unsigned ix, iy, iz;
  double x, y, z;
  const Vector& xmin = outerBoundary.xmin();
  const Vector& xmax = outerBoundary.xmax();
  const double dx = (xmax - xmin).maxElement()/n1d;
  VERIFY(dx > 0.0);
  const unsigned nx = std::max(1U, unsigned((xmax.x() - xmin.x())/dx + 0.5));
  const unsigned ny = std::max(1U, unsigned((xmax.y() - xmin.y())/dx + 0.5));
  const unsigned nz = std::max(1U, unsigned((xmax.z() - xmin.z())/dx + 0.5));
  CHECK(nx > 0 and nx <= n1d);
  CHECK(ny > 0 and ny <= n1d);
  CHECK(nz > 0 and nz <= n1d);
  result.reserve(nx*ny*nz);

  if (outerBoundary.convex()) {

    for (iz = 0; iz != nz; ++iz) {
      z = xmin.z() + (iz + 0.5)*dx;
      for (iy = 0; iy != ny; ++iy) {
        y = xmin.y() + (iy + 0.5)*dx;
        for (ix = 0; ix != nx; ++ix) {
          x = xmin.x() + (ix + 0.5)*dx;
          Vector pos(x, y, z);
          if (outerBoundary.contains(pos)) result.push_back(pos);
        }
      }
    }

  } else {

    // We'll try using the much faster convex contains method as a pre-filter
    // before using the expensive general contains.
    const FacetedVolume convexOuterBoundary(outerBoundary.vertices());
    for (iz = 0; iz != nz; ++iz) {
      z = xmin.z() + (iz + 0.5)*dx;
      for (iy = 0; iy != ny; ++iy) {
        y = xmin.y() + (iy + 0.5)*dx;
        for (ix = 0; ix != nx; ++ix) {
          x = xmin.x() + (ix + 0.5)*dx;
          Vector pos(x, y, z);
          if (convexOuterBoundary.convexContains(pos) and outerBoundary.contains(pos)) result.push_back(pos);
        }
      }
    }

  }
  return result;
}

//------------------------------------------------------------------------------
// Fill between an inner and outer boundary.
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& innerBoundary,
                  const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  vector<Vector> result;
  unsigned ix, iy, iz;
  double x, y, z;
  const Vector& xmin = outerBoundary.xmin();
  const Vector& xmax = outerBoundary.xmax();
  const double dx = (xmax - xmin).maxElement()/n1d;
  VERIFY(dx > 0.0);
  const unsigned nx = std::max(1U, unsigned((xmax.x() - xmin.x())/dx + 0.5));
  const unsigned ny = std::max(1U, unsigned((xmax.y() - xmin.y())/dx + 0.5));
  const unsigned nz = std::max(1U, unsigned((xmax.z() - xmin.z())/dx + 0.5));
  CHECK(nx > 0 and nx <= n1d);
  CHECK(ny > 0 and ny <= n1d);
  CHECK(nz > 0 and nz <= n1d);
  result.reserve(nx*ny*nz);

  const bool innerConvex = innerBoundary.convex();
  const bool outerConvex = outerBoundary.convex();
  const FacetedVolume convexInnerBoundary = innerConvex ? innerBoundary : FacetedVolume(innerBoundary.vertices());
  const FacetedVolume convexOuterBoundary = outerConvex ? outerBoundary : FacetedVolume(outerBoundary.vertices());

  bool outerTest, innerTest;
  for (iz = 0; iz != nz; ++iz) {
    z = xmin.z() + (iz + 0.5)*dx;
    for (iy = 0; iy != ny; ++iy) {
      y = xmin.y() + (iy + 0.5)*dx;
      for (ix = 0; ix != nx; ++ix) {
        x = xmin.x() + (ix + 0.5)*dx;
        Vector pos(x, y, z);
        outerTest = convexOuterBoundary.contains(pos);
        if (outerTest and not outerConvex) outerTest = outerBoundary.contains(pos);
        if (outerTest) {
          innerTest = convexInnerBoundary.contains(pos);
          if (not innerTest) {
            result.push_back(pos);
          } else {
            if (not innerConvex) innerTest = innerBoundary.contains(pos);
            if (not innerTest) result.push_back(pos);
          }
        }
      }
    }
  }

  return result;
}

}
