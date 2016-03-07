//------------------------------------------------------------------------------
// pointOnPolyhedron
//
// Test if a given point is on the boundary of a polyhedron.
// You should do the bounding box check before calling this method for 
// efficiency!  We do not repeat that check here since this is meant to be 
// called by code that has already done that check.
//------------------------------------------------------------------------------
#ifndef __Spheral_pointOnPolyhedron__
#define __Spheral_pointOnPolyhedron__

#include "Geometry/Dimension.hh"

namespace Spheral {

bool pointOnPolyhedron(const Dim<3>::Vector& p,
                       const Dim<3>::FacetedVolume& polyhedron,
                       const double tol = 1.0e-10);

}

#endif
