text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeGenerators/compactFacetedVolumes.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template
unsigned compactFacetedVolumes<Dim<%(ndim)s>>(std::vector<Dim<%(ndim)s>::FacetedVolume>& shapes,
                                              std::vector<Dim<%(ndim)s>::Vector>& centers,
                                              std::vector<int>& flags,
                                              const Dim<%(ndim)s>::FacetedVolume& surface,
                                              const double depthmax,
                                              const unsigned surfaceIterations,
                                              const unsigned maxIterations,
                                              const double dispfrac,
                                              const double maxoverlapfrac);
}
"""
