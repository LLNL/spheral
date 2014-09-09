//---------------------------------Spheral++----------------------------------//
// refinePolyhedron
//
// Refine/smooth a polyhedron.
//
// Created by JMO, Tue Mar 18 15:23:04 PDT 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_refinePolyhedron__
#define __Spheral_refinePolyhedron__

#include "Geometry/GeomPolyhedron.hh"

namespace Spheral {
  GeomPolyhedron refinePolyhedron(const GeomPolyhedron& poly0,
                                  const unsigned numLevels);
}

#endif

