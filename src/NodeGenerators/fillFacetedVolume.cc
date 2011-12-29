#include "fillFacetedVolume.hh"

namespace Spheral {

using namespace std;
using std::abs;
using std::min;
using std::max;

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

  for (iz = 0; iz != nz; ++iz) {
    z = xmin.z() + (iz + 0.5)*dx;
    for (iy = 0; iy != ny; ++iy) {
      y = xmin.y() + (iy + 0.5)*dx;
      for (ix = 0; ix != nx; ++ix) {
        x = xmin.x() + (ix + 0.5)*dx;
        if (outerBoundary.contains(Vector(x, y, z))) result.push_back(Vector(x, y, z));
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

  for (iz = 0; iz != nz; ++iz) {
    z = xmin.z() + (iz + 0.5)*dx;
    for (iy = 0; iy != ny; ++iy) {
      y = xmin.y() + (iy + 0.5)*dx;
      for (ix = 0; ix != nx; ++ix) {
        x = xmin.x() + (ix + 0.5)*dx;
        Vector pos(x, y, z);
        if ((not innerBoundary.contains(pos)) and outerBoundary.contains(pos)) result.push_back(pos);
      }
    }
  }

  return result;
}

}
