#include "fillFacetedVolume.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {
//------------------------------------------------------------------------------
// Compute an individual HCP position.
//------------------------------------------------------------------------------
inline
void HCPposition(const unsigned i,
                 const unsigned nxy,
                 const unsigned nx,
                 const unsigned ny,
                 const double dx,
                 const double dy,
                 const double dz,
                 const Dim<3>::Vector& xmin,
                 const Dim<3>::Vector& xmax,
                 Dim<3>::Vector& result) {
  const auto ix = i % nx;
  const auto iy = (i / nx) % ny;
  const auto iz = i / nxy;
  result[0] = xmin[0] + (ix + 0.5*((iy % 2) + (iz % 2)))*dx;
  result[1] = xmin[1] + (iy + 0.5*(iz % 2))*dy;
  result[2] = xmin[2] + (iz + 0.5)*dz;
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
  unsigned ix, iy, iz;
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
#pragma omp parallel
  {
    Vector pos;
    vector<Vector> result_local;
    result_local.reserve(ndomain);
#pragma omp for
    for (auto i = imin; i < imax; ++i) {
      HCPposition(i, nxy, nx, ny, dx, dx, dx, xmin, xmax, pos);    
      if (outerBoundary.contains(pos, false)) result_local.push_back(pos);
    }
#pragma omp critical
    {
      result.insert(result.end(), result_local.begin(), result_local.end());
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Fill an outer bounding volume where dx is user-defined.
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume2(const Dim<3>::FacetedVolume& outerBoundary,
                   const double   dx,
                   const unsigned domain,
                   const unsigned numDomains) {
  VERIFY(dx > 0.0);
  VERIFY(numDomains >= 1);
  VERIFY(domain < numDomains);

  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  vector<Vector> result;
  unsigned ix, iy, iz;
  const Vector& xmin = outerBoundary.xmin();
  const Vector& xmax = outerBoundary.xmax();
  const unsigned nx = std::max(1U, unsigned((xmax.x() - xmin.x())/dx + 0.5));
  const unsigned ny = std::max(1U, unsigned((xmax.y() - xmin.y())/dx + 0.5));
  const unsigned nz = std::max(1U, unsigned((xmax.z() - xmin.z())/dx + 0.5));
  CHECK(nx > 0);
  CHECK(ny > 0);
  CHECK(nz > 0);

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
#pragma omp parallel
  {
    Vector pos;
    vector<Vector> result_local;
    result_local.reserve(ndomain);
#pragma omp for
    for (auto i = imin; i < imax; ++i) {
      HCPposition(i, nxy, nx, ny, dx, dx, dx, xmin, xmax, pos);    
      if (outerBoundary.contains(pos, false)) result.push_back(pos);
    }
#pragma omp critical
    {
      result.insert(result.end(), result_local.begin(), result_local.end());
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
                  const unsigned n1d,
                  const unsigned domain,
                  const unsigned numDomains) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;

  vector<Vector> result;
  unsigned ix, iy, iz;
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

#pragma omp parallel
  {
    Vector pos;
    vector<Vector> result_local;
    result_local.reserve(ndomain);
#pragma omp for
    for (auto i = imin; i < imax; ++i) {
      HCPposition(i, nxy, nx, ny, dx, dx, dx, xmin, xmax, pos);    
      if (outerBoundary.contains(pos, false) and not innerBoundary.contains(pos, false)) {
        result.push_back(pos);
      }
    }
#pragma omp critical
    {
      result.insert(result.end(), result_local.begin(), result_local.end());
    }
  }
  return result;
}

}
