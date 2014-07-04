//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on convex hulls.
//------------------------------------------------------------------------------
#include "computeHullVolumes.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

using namespace std;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using NodeSpace::NodeList;

template<typename Dimension>
void
computeHullVolumes(const ConnectivityMap<Dimension>& connectivityMap,
                   const FieldList<Dimension, typename Dimension::Vector>& position,
                   FieldList<Dimension, typename Dimension::Scalar>& volume) {

  // Pre-conditions.
  const size_t numNodeLists = volume.size();
  REQUIRE(position.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Zero out the result.
  volume = 0.0;

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);

      // Collect the positions of all neighbors.
      vector<Vector> positionsInv(1, Vector::zero);
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;
          const Vector& rj = position(nodeListj, j);
          const Vector rji = rj - ri,
                       rjiHat = rji.unitVector();
          positionsInv.push_back(1.0/sqrt(rji.magnitude2() + 1.0e-30) * rjiHat);
        }
      }

      // Build the hull of the inverse.
      const FacetedVolume hullInv(positionsInv);

      // Use the vertices selected by the inverse hull to construct the
      // volume of the node.
      vector<Vector> positions;
      const vector<Vector>& vertsInv = hullInv.vertices();
      for (typename std::vector<Vector>::const_iterator itr = vertsInv.begin();
           itr != vertsInv.end();
           ++itr) {
        positions.push_back(1.0/sqrt(itr->magnitude2() + 1.0e-30) * itr->unitVector());
      }
      const FacetedVolume hulli(positions, hullInv.facetVertices());

      // And now we have the volume.
      volume(nodeListi, i) = hulli.volume();
    }
  }
}

}

