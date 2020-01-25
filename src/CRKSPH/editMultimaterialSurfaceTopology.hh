//------------------------------------------------------------------------------
// Look for any points that touch a surface (multi-material or void).
// For such points:
//   - Remove any non-surface multi-material overlap.
//   - If not a surface point, flag this point as touching a surface point with
//     surfacePoint=-1.
//------------------------------------------------------------------------------
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

template<typename Dimension>
void
editMultimaterialSurfaceTopology(FieldList<Dimension, int>& surfacePoint,
                                 ConnectivityMap<Dimension>& connectivityMap);

}
