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
// Fill an outer bounding volume (specify x number of points).
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d,
                  const unsigned domain,
                  const unsigned numDomains) {
  VERIFY(n1d > 0);
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  const auto& xmin = outerBoundary.xmin();
  const auto& xmax = outerBoundary.xmax();
  const auto  dx = (xmax - xmin).maxElement()/n1d;
  return fillFacetedVolume10(outerBoundary, FacetedVolume(), dx, domain, numDomains);
}

//------------------------------------------------------------------------------
// Fill an outer bounding volume (dx specified).
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume2(const Dim<3>::FacetedVolume& outerBoundary,
                   const double dx,
                   const unsigned domain,
                   const unsigned numDomains) {
  VERIFY(dx > 0.0);
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  return fillFacetedVolume10(outerBoundary, FacetedVolume(), dx, domain, numDomains);
}

//------------------------------------------------------------------------------
// Fill between an inner and outer boundary (specify x number of points).
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume3(const Dim<3>::FacetedVolume& innerBoundary,
                   const Dim<3>::FacetedVolume& outerBoundary,
                   const unsigned n1d,
                   const unsigned domain,
                   const unsigned numDomains) {
  VERIFY(n1d > 0);
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  const auto& xmin = outerBoundary.xmin();
  const auto& xmax = outerBoundary.xmax();
  const auto  dx = (xmax - xmin).maxElement()/n1d;
  return fillFacetedVolume10(outerBoundary, FacetedVolume(), dx, domain, numDomains);
}

//------------------------------------------------------------------------------
// Fill between an inner and outer boundary (dx specified).
//------------------------------------------------------------------------------
vector<Dim<3>::Vector>
fillFacetedVolume10(const Dim<3>::FacetedVolume& outerBoundary,
                    const Dim<3>::FacetedVolume& innerBoundary,
                    const double dx,
                    const unsigned domain,
                    const unsigned numDomains) {
  typedef Dim<3>::Vector Vector;
  VERIFY(dx > 0.0);
  VERIFY(numDomains >= 1);
  VERIFY(domain < numDomains);

  const bool useInner = not (innerBoundary.facets().empty());
  vector<Vector> result;

  // Find numbers of points
  const auto& xmin = outerBoundary.xmin();
  const auto& xmax = outerBoundary.xmax();
  const auto nx = std::max(1U, unsigned((xmax.x() - xmin.x())/dx + 0.5));
  const auto ny = std::max(1U, unsigned((xmax.y() - xmin.y())/dx + 0.5));
  const auto nz = std::max(1U, unsigned((xmax.z() - xmin.z())/dx + 0.5));
  const Vector delta((xmax[0] - xmin[0])/nx,
                     (xmax[1] - xmin[1])/ny,
                     (xmax[2] - xmin[2])/nz);

  // Pick the longest direction to construct our ray.
  const auto length = 2.0*(xmax - xmin).magnitude();
  unsigned stride_index, iindex, jindex, nplanerays1, nplanerays2;
  Vector ray1, ray2;
  if (nx >= ny and nx >= nz) {
    nplanerays1 = ny;
    nplanerays2 = nz;
    stride_index = 0;
    iindex = 1;
    jindex = 2;
  } else if (ny >= nx and ny >= nz) {
    nplanerays1 = nx;
    nplanerays2 = nz;
    iindex = 0;
    stride_index = 1;
    jindex = 2;
  } else {
    nplanerays1 = nx;
    nplanerays2 = ny;
    iindex = 0;
    jindex = 1;
    stride_index = 2;
  }
  ray1[stride_index] = xmin[stride_index] - 0.1*length;
  ray2[stride_index] = xmax[stride_index] + 0.1*length;
  const auto nplanerays = nplanerays1*nplanerays2;

  // How should we carve up the ray-casting work between MPI domains?
  const unsigned ndomain0 = nplanerays/numDomains;
  const unsigned remainder = nplanerays - ndomain0*numDomains;
  CHECK(remainder < numDomains);
  const unsigned ndomain = nplanerays/numDomains + (domain < remainder ? 1 : 0);
  const unsigned imin = domain*ndomain0 + min(domain, remainder);
  const unsigned imax = imin + ndomain;
  CHECK(domain < numDomains - 1 or imax == nplanerays);

  // Walk the projection plane axes.
#pragma omp parallel
  {
    vector<unsigned> facetIDs;
    vector<Vector> intersections;
    vector<Vector> result_thread;
#pragma omp for firstprivate (ray1, ray2)
    for (auto nray = imin; nray < imax; ++nray) {
      auto iray = nray % nplanerays1;
      auto jray = nray / nplanerays1;
      ray1[iindex] = xmin[iindex] + (iray + 0.5)*delta[iindex];
      ray1[jindex] = xmin[jindex] + (jray + 0.5)*delta[jindex];
      ray2[iindex] = ray1[iindex];
      ray2[jindex] = ray1[jindex];

      // Project through the volume and find intersection points (should be in pairs).
      outerBoundary.intersect(ray1, ray2, facetIDs, intersections);
      const auto nintersect = intersections.size();
      CHECK(facetIDs.size() == nintersect);
      // CHECK(nintersect % 2 == 0);

      // Now points between pairs of intersections should be interior to the surface.
      for (auto i = 0; i < nintersect;) {
        if (i + 1 < nintersect) {
          const auto& x1 = intersections[i];
          const auto& x2 = intersections[i + 1];

          // Because of degeneracies, it's possible to get isolated intersection points
          // rather than a pair bracketing interior space.  We have to check...
          if (outerBoundary.contains(0.5*(x1 + x2), false)) {
            const auto  ni = max(1, int((x2[stride_index] - x1[stride_index])/dx + 0.5));
            const auto  dstep = (x2 - x1)/ni;
            for (auto k = 0; k < ni; ++k) result_thread.push_back(x1 + (k + 0.5)*dstep);
          }
          i += 2;    // This pair was interior, so skip to the next possible start
        } else {
          i += 1;    // Looking for the next pair bounding interior
        }
      }
    }
#pragma omp critical
    {
      result.insert(result.end(), result_thread.begin(), result_thread.end());
    }
  }
  return result;
}

}
