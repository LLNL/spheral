//---------------------------------Spheral++------------------------------------
// Given our integer FieldList of surface flags, spread that info to the 
// neighboring points.
//------------------------------------------------------------------------------
#include "flagSurfaceNeighbors.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NeighborSpace::ConnectivityMap;

template<typename Dimension>
void
flagSurfaceNeighbors(FieldSpace::FieldList<Dimension, int>& surfacePoint,
                     const NeighborSpace::ConnectivityMap<Dimension >& connectivityMap) {

  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::SymTensor SymTensor;
  typedef Dim<1>::FacetedVolume FacetedVolume;

  const unsigned numNodeLists = surfacePoint.size();

  // Walk the NodeLists.
  for (unsigned nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const unsigned n = surfacePoint[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {

      // If this point is as yet unflagged, check for any surface neighbors.
      if (surfacePoint(nodeListi, i) == 0) {
        const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          for (vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
               jItr != fullConnectivity[nodeListj].end();
               ++jItr) {
            const unsigned j = *jItr;
            if (surfacePoint(nodeListj, j) > 0) {
              surfacePoint(nodeListi, i) = -1;
              continue;
            }
          }
        }
      }
    }
  }
}

}
}
