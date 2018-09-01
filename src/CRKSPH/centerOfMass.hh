//---------------------------------Spheral++------------------------------------
// Compute the center of mass of a FacetedVolume assuming a linear mass density
// field.
//------------------------------------------------------------------------------
#ifndef __Spheral__centerOfMass__
#define __Spheral__centerOfMass__

#include "Geometry/Dimension.hh"

namespace Spheral {
  Dim<1>::Vector centerOfMass(const Dim<1>::FacetedVolume& polyvol,
                              const Dim<1>::Vector& gradRhoi);
  Dim<2>::Vector centerOfMass(const Dim<2>::FacetedVolume& polyvol,
                              const Dim<2>::Vector& gradRhoi);
  Dim<3>::Vector centerOfMass(const Dim<3>::FacetedVolume& polyvol,
                              const Dim<3>::Vector& gradRhoi);
}

#endif
