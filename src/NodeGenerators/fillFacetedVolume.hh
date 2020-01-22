//------------------------------------------------------------------------------
// Return a set of positions filling a faceted volume, optionally with an
// inner boundary as well.
//------------------------------------------------------------------------------
#ifndef __Spheral_fillFacetedVolume__
#define __Spheral_fillFacetedVolume__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

// Fill a bounding volume.
std::vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d,
                  const unsigned domain,
                  const unsigned numDomains);

// Fill a bounding volume where dx is user-defined.
std::vector<Dim<3>::Vector>
fillFacetedVolume2(const Dim<3>::FacetedVolume& outerBoundary,
                   const double   dx,
                   const unsigned domain,
                   const unsigned numDomains);

// Fill between inner and outer bounding volumes.
std::vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& innerBoundary,
                  const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d,
                  const unsigned domain,
                  const unsigned numDomains);

}

#endif
