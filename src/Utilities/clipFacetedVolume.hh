//------------------------------------------------------------------------------
// A collection of Spheral wrappers for using PolyClipper methods.
//------------------------------------------------------------------------------
#ifndef __Spheral_clipFacetedVolume__
#define __Spheral_clipFacetedVolume__

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Clip a polygon/polyhedron by a series of planes
//------------------------------------------------------------------------------
Dim<2>::FacetedVolume clipFacetedVolume(const Dim<2>::FacetedVolume& poly,
                                        const std::vector<GeomPlane<Dim<2> > >& planes);

Dim<3>::FacetedVolume clipFacetedVolume(const Dim<3>::FacetedVolume& poly,
                                        const std::vector<GeomPlane<Dim<3> > >& planes);

//------------------------------------------------------------------------------
// Return the volume of the clipped region.
//------------------------------------------------------------------------------
double clippedVolume(const Dim<2>::FacetedVolume& poly,
                     const std::vector<GeomPlane<Dim<2> > >& planes);

double clippedVolume(const Dim<3>::FacetedVolume& poly,
                     const std::vector<GeomPlane<Dim<3> > >& planes);

}

#endif
