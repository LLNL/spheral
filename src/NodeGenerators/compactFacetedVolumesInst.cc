//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeGenerators/compactFacetedVolumes.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template
unsigned compactFacetedVolumes<Dim<1>>(std::vector<Dim<1>::FacetedVolume>& shapes,
                                              std::vector<Dim<1>::Vector>& centers,
                                              std::vector<int>& flags,
                                              const Dim<1>::FacetedVolume& surface,
                                              const double depthmax,
                                              const unsigned surfaceIterations,
                                              const unsigned maxIterations,
                                              const double dispfrac,
                                              const double maxoverlapfrac);
#endif

#if defined(SPHERAL_ENABLE_2D)
template
unsigned compactFacetedVolumes<Dim<2>>(std::vector<Dim<2>::FacetedVolume>& shapes,
                                              std::vector<Dim<2>::Vector>& centers,
                                              std::vector<int>& flags,
                                              const Dim<2>::FacetedVolume& surface,
                                              const double depthmax,
                                              const unsigned surfaceIterations,
                                              const unsigned maxIterations,
                                              const double dispfrac,
                                              const double maxoverlapfrac);
#endif

#if defined(SPHERAL_ENABLE_3D)
template
unsigned compactFacetedVolumes<Dim<3>>(std::vector<Dim<3>::FacetedVolume>& shapes,
                                              std::vector<Dim<3>::Vector>& centers,
                                              std::vector<int>& flags,
                                              const Dim<3>::FacetedVolume& surface,
                                              const double depthmax,
                                              const unsigned surfaceIterations,
                                              const unsigned maxIterations,
                                              const double dispfrac,
                                              const double maxoverlapfrac);
#endif
}