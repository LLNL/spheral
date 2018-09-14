text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeGenerators/chooseRandomNonoverlappingCenter.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template
unsigned chooseRandomNonoverlappingCenter<Dim<%(ndim)s>>(Dim<%(ndim)s>::Vector& result,
                                 const Dim<%(ndim)s>::FacetedVolume& shape,
                                 const Dim<%(ndim)s>::FacetedVolume& boundary,
                                 const std::vector<Dim<%(ndim)s>::FacetedVolume>& existingShapes,
                                 const unsigned maxTries);
}
"""
