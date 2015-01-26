#include "fillFacetedVolume.hh"

namespace Spheral {

using namespace std;

namespace {
//------------------------------------------------------------------------------
// Compute an individual HCP position.
//------------------------------------------------------------------------------
inline
Dim<3>::Vector HCPposition(const unsigned i,
                           const unsigned nx,
                           const unsigned ny,
                           const unsigned nz,
                           const double dx,
                           const double dy,
                           const double dz,
                           const Dim<3>::Vector& xmin,
                           const Dim<3>::Vector& xmax) {
  const unsigned nxy = nx*ny;
  const unsigned ix = i % nx;
  const unsigned iy = (i / nx) % ny;
  const unsigned iz = i / nxy;
  return Dim<3>::Vector(xmin.x() + (ix + 0.5*((iy % 2) + (iz % 2)))*dx,
                        xmin.y() + (iy + 0.5*(iz % 2))*dy,
                        xmin.z() + (iz + 0.5)*dz);
}
}

//------------------------------------------------------------------------------
// Fill an outer bounding volume.
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d,
                  const unsigned domain,
                  const unsigned numDomains) {
  VERIFY(n1d > 0);
  VERIFY(numDomains >= 1);
  VERIFY(domain < numDomains);

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  vector<Vector> result;
  unsigned i, ix, iy, iz;
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

  // Figure out the range of global IDs for this domain.
  const unsigned nxy = nx*ny;
  const unsigned nxyz = nx*ny*nz;
  const unsigned ndomain0 = nxyz/numDomains;
  const unsigned remainder = nxyz - ndomain0*numDomains;
  CHECK(remainder < numDomains);
  const unsigned ndomain = nxyz/numDomains + (domain < remainder ? 1 : 0);
  const unsigned imin = domain*ndomain0 + min(domain, remainder);
  const unsigned imax = imin + ndomain;
  CHECK(domain < numDomains - 1 or imax == nxyz);
  result.reserve(ndomain);

  // Use an HCP lattice to do the filling.
  for (i = imin; i != imax; ++i) {
    const Vector pos = HCPposition(i, nx, ny, nz, dx, dx, dx, xmin, xmax);    
    if (outerBoundary.contains(pos)) result.push_back(pos);
  }
  return result;
}

//------------------------------------------------------------------------------
// Fill between an inner and outer boundary.
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& innerBoundary,
                  const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d,
                  const unsigned domain,
                  const unsigned numDomains) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  vector<Vector> result;
  unsigned i, ix, iy, iz;
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

  // Figure out the range of global IDs for this domain.
  const unsigned nxy = nx*ny;
  const unsigned nxyz = nx*ny*nz;
  const unsigned ndomain0 = nxyz/numDomains;
  const unsigned remainder = nxyz - ndomain0*numDomains;
  CHECK(remainder < numDomains);
  const unsigned ndomain = nxyz/numDomains + (domain < remainder ? 1 : 0);
  const unsigned imin = domain*ndomain0 + min(domain, remainder);
  const unsigned imax = imin + ndomain;
  CHECK(domain < numDomains - 1 or imax == nxyz);
  result.reserve(ndomain);

  for (i = imin; i != imax; ++i) {
    const Vector pos = HCPposition(i, nx, ny, nz, dx, dx, dx, xmin, xmax);    
    if (outerBoundary.contains(pos) and not innerBoundary.contains(pos)) {
        result.push_back(pos);
    }
  }

  return result;
}

}
