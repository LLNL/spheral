//---------------------------------Spheral++----------------------------------//
// PolyClipperUtilities
//
// Methods for talking with PolyClipper:
//       https://github.com/LLNL/PolyClipper
//
// Created by JMO, Thu Feb  4 16:24:44 PST 2021
//----------------------------------------------------------------------------//
#ifndef __PolyClipperUtilities__
#define __PolyClipperUtilities__

#include "GeomVector.hh"
#include "GeomPolygon.hh"
#include "GeomPolyhedron.hh"
#include "PolyClipper/polyclipper.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Polygons
//------------------------------------------------------------------------------
void convertToPolyClipper(PolyClipper::Polygon& polygon,
                          const Dim<2>::FacetedVolume& Spheral_polygon);

std::vector<std::set<int>> convertFromPolyClipper(Dim<2>::FacetedVolume& Spheral_polygon,
                                                  const PolyClipper::Polygon& polygon);

//------------------------------------------------------------------------------
// Polyhedra
//------------------------------------------------------------------------------
void convertToPolyClipper(PolyClipper::Polyhedron& polyhedron,
                          const Dim<2>::FacetedVolume& Spheral_polyhedron);

std::vector<std::set<int>> convertFromPolyClipper(Dim<2>::FacetedVolume& Spheral_polyhedron,
                                                  const PolyClipper::Polyhedron& polyhedron);

}

#endif
