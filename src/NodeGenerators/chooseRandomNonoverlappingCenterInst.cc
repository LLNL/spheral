//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeGenerators/chooseRandomNonoverlappingCenter.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template
unsigned chooseRandomNonoverlappingCenter<Dim<1>>(Dim<1>::Vector& result,
                                 const Dim<1>::FacetedVolume& shape,
                                 const Dim<1>::FacetedVolume& boundary,
                                 const std::vector<Dim<1>::FacetedVolume>& existingShapes,
                                 const unsigned maxTries);
#endif

#if defined(SPHERAL_ENABLE_2D)
template
unsigned chooseRandomNonoverlappingCenter<Dim<2>>(Dim<2>::Vector& result,
                                 const Dim<2>::FacetedVolume& shape,
                                 const Dim<2>::FacetedVolume& boundary,
                                 const std::vector<Dim<2>::FacetedVolume>& existingShapes,
                                 const unsigned maxTries);
#endif

#if defined(SPHERAL_ENABLE_3D)
template
unsigned chooseRandomNonoverlappingCenter<Dim<3>>(Dim<3>::Vector& result,
                                 const Dim<3>::FacetedVolume& shape,
                                 const Dim<3>::FacetedVolume& boundary,
                                 const std::vector<Dim<3>::FacetedVolume>& existingShapes,
                                 const unsigned maxTries);
#endif
}