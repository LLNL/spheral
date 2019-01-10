//---------------------------------Spheral++----------------------------------//
// copy2polytope
//
// Helper method to copy a set of Spheral polyhedra to a polytope tessellation.
//
// Created by JMO, Wed Jan  9 16:25:28 PST 2019
//----------------------------------------------------------------------------//
#ifndef __Spheral_copy2polytope_hh__
#define __Spheral_copy2polytope_hh__

#include "Field/FieldList.hh"
#include "polytope/Tessellation.hh"

namespace Spheral {

template<typename Dimension, int nDim>
void copy2polytope(const FieldList<Dimension, typename Dimension::FacetedVolume>& cells,
                   polytope::Tessellation<nDim, double>& mesh);
                   
}

#endif

