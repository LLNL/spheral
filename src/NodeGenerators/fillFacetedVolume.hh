//------------------------------------------------------------------------------
// Return a set of positions filling a faceted volume, optionally with an
// inner boundary as well.
//------------------------------------------------------------------------------
#ifndef __Spheral_fillFacetedVolume__
#define __Spheral_fillFacetedVolume__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

// Fill an outer bounding volume (specify x number of points).
std::vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d,
                  const unsigned domain,
                  const unsigned numDomains);

// Fill an outer bounding volume (dx specified).
std::vector<Dim<3>::Vector>
fillFacetedVolume2(const Dim<3>::FacetedVolume& outerBoundary,
                   const double dx,
                   const unsigned domain,
                   const unsigned numDomains);

// Fill between an inner and outer boundary (specify x number of points).
std::vector<Dim<3>::Vector>
fillFacetedVolume3(const Dim<3>::FacetedVolume& innerBoundary,
                   const Dim<3>::FacetedVolume& outerBoundary,
                   const unsigned n1d,
                   const unsigned domain,
                   const unsigned numDomains);

// Fill between an inner and outer boundary (dx specified).
std::vector<Dim<3>::Vector>
fillFacetedVolume10(const Dim<3>::FacetedVolume& innerBoundary,
                    const Dim<3>::FacetedVolume& outerBoundary,
                    const double dx,
                    const unsigned domain,
                    const unsigned numDomains);

}

#endif
