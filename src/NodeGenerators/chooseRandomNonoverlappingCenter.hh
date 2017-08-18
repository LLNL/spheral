#ifndef __Spheral_chooseRandomNonoverlappingCenter__
#define __Spheral_chooseRandomNonoverlappingCenter__

#include <vector>

//------------------------------------------------------------------------------
// Push FacetedVolume shapes together inside a surface, but excluding mutual
// overlap.
//------------------------------------------------------------------------------
namespace Spheral {

template<typename Dimension>
unsigned
chooseRandomNonoverlappingCenter(typename Dimension::Vector& result,
                                 const typename Dimension::FacetedVolume& trialShape,
                                 const typename Dimension::FacetedVolume& boundary,
                                 const std::vector<typename Dimension::FacetedVolume>& existingShapes,
                                 const unsigned maxTries);

}

#endif
