#ifndef __Spheral_compactFacetedVolumes__
#define __Spheral_compactFacetedVolumes__

#include <vector>

//------------------------------------------------------------------------------
// Push FacetedVolume shapes together inside a surface, but excluding mutual
// overlap.
//------------------------------------------------------------------------------
namespace Spheral {
template<typename Dimension>
unsigned compactFacetedVolumes(std::vector<typename Dimension::FacetedVolume>& shapes,
                               std::vector<typename Dimension::Vector>& centers,
                               std::vector<int>& flags,
                               const typename Dimension::FacetedVolume& surface,
                               const double depthmax,
                               const unsigned surfaceIterations,
                               const unsigned maxIterations,
                               const double dispfrac,
                               const double maxoverlapfrac);
}

#endif
