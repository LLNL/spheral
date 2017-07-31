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
// Convert a Spheral polygon/polyhedron to a r2d_poly/r3d_poly to
//------------------------------------------------------------------------------
void polygon_to_r2d_poly(const Dim<2>::FacetedVolume& poly, r2d_poly& result);
void polyhedron_to_r3d_poly(const Dim<3>::FacetedVolume& poly, r3d_poly& result);

//------------------------------------------------------------------------------
// Convert a r2d_poly/r3d_poly to a Spheral polygon/polyhedron.
//------------------------------------------------------------------------------
void r2d_poly_to_polygon(const r2d_poly& celli,
                         const double tol,
                         Dim<2>::FacetedVolume& result);
void r3d_poly_to_polyhedron(const r3d_poly& celli,
                            const double tol,
                            Dim<3>::FacetedVolume& result);

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
