//---------------------------------Spheral++----------------------------------//
// coarsenBinnedValues
//
// Given a lattice of Values, return the multi-level result of progressively
// coarsening the distribution a requeted number of levels.  The passed 
// in vector<vector<Value> > should have the finest values as the last element,
// and be sized the required number of levels you want to coarsen.  This method
// simply changes that vector<vector> in place.
//
// Created by JMO, Fri Feb 19 09:38:29 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_coarsenBinnedValues__
#define __Spheral_coarsenBinnedValues__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// 1-D.
//------------------------------------------------------------------------------
template<typename Value>
void
coarsenBinnedValues(std::vector<std::vector<Value> >& values,
                    const unsigned nxFine);

//------------------------------------------------------------------------------
// 2-D.
//------------------------------------------------------------------------------
template<typename Value>
void
coarsenBinnedValues(std::vector<std::vector<Value> >& values,
                    const unsigned nxFine,
                    const unsigned nyFine);

//------------------------------------------------------------------------------
// 3-D.
//------------------------------------------------------------------------------
template<typename Value>
void
coarsenBinnedValues(std::vector<std::vector<Value> >& values,
                    const unsigned nxFine,
                    const unsigned nyFine,
                    const unsigned nzFine);

}

#endif
