#include "fillFacetedVolume.hh"

namespace Spheral {

using namespace std;

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

  if (outerBoundary.convex()) {

    for (i = imin; i != imax; ++i) {
      ix = i % nx;
      iy = (i / nx) % ny;
      iz = i / nxy;
      const Vector pos(xmin.x() + (ix + 0.5)*dx,
                       xmin.y() + (iy + 0.5)*dx,
                       xmin.z() + (iz + 0.5)*dx);
      if (outerBoundary.contains(pos)) result.push_back(pos);
    }

  } else {

    // We'll try using the much faster convex contains method as a pre-filter
    // before using the expensive general contains.
    const FacetedVolume convexOuterBoundary(outerBoundary.vertices());
    for (i = imin; i != imax; ++i) {
      ix = i % nx;
      iy = (i / nx) % ny;
      iz = i / nxy;
      const Vector pos(xmin.x() + (ix + 0.5)*dx,
                       xmin.y() + (iy + 0.5)*dx,
                       xmin.z() + (iz + 0.5)*dx);
      if (convexOuterBoundary.convexContains(pos) and outerBoundary.contains(pos)) result.push_back(pos);
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

  const bool innerConvex = innerBoundary.convex();
  const bool outerConvex = outerBoundary.convex();
  const FacetedVolume convexInnerBoundary = innerConvex ? innerBoundary : FacetedVolume(innerBoundary.vertices());
  const FacetedVolume convexOuterBoundary = outerConvex ? outerBoundary : FacetedVolume(outerBoundary.vertices());

  bool outerTest, innerTest;
  for (i = imin; i != imax; ++i) {
    ix = i % nx;
    iy = (i / nx) % ny;
    iz = i / nxy;
    const Vector pos(xmin.x() + (ix + 0.5)*dx,
                     xmin.y() + (iy + 0.5)*dx,
                     xmin.z() + (iz + 0.5)*dx);
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

  return result;
}

}
