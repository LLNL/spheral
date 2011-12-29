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
                  const unsigned n1d);

// Fill between inner and outer bounding volumes.
std::vector<Dim<3>::Vector>
fillFacetedVolume(const Dim<3>::FacetedVolume& innerBoundary,
                  const Dim<3>::FacetedVolume& outerBoundary,
                  const unsigned n1d);

}

#endif
