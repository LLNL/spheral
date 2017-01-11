//------------------------------------------------------------------------------
// A collection of Spheral wrappers for using R2D/R3D methods.
//------------------------------------------------------------------------------
#ifndef __Spheral_r3d_utils__
#define __Spheral_r3d_utils__

extern "C" {
#include "r3d/r2d.h"
#include "r3d/r3d.h"
}

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Convert an r2d_poly/r3d_poly to a Spheral polygon/polyhedron.
//------------------------------------------------------------------------------
Dim<2>::FacetedVolume r2d_poly_to_polygon(const r2d_poly& celli,
                                          const double tol);

Dim<3>::FacetedVolume r3d_poly_to_polyhedron(const r3d_poly& celli,
                                             const double tol);

//------------------------------------------------------------------------------
// Clip a polygon/polyhedron by a series of planes
//------------------------------------------------------------------------------
Dim<2>::FacetedVolume clipFacetedVolume(const Dim<2>::FacetedVolume& poly,
                                        const std::vector<GeomPlane<Dim<2> > >& planes);

Dim<3>::FacetedVolume clipFacetedVolume(const Dim<3>::FacetedVolume& poly,
                                        const std::vector<GeomPlane<Dim<3> > >& planes);

}

#endif
