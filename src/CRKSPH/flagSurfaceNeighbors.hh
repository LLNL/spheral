//---------------------------------Spheral++------------------------------------
// Given our integer FieldList of surface flags, spread that info to the 
// neighboring points.
//------------------------------------------------------------------------------
#ifndef __Spheral__flagSurfaceNeighbors__
#define __Spheral__flagSurfaceNeighbors__

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {
template<typename Dimension> 
void flagSurfaceNeighbors(FieldList<Dimension, int>& surfacePoint,
                          const ConnectivityMap<Dimension>& connectivityMap);
}

#endif
