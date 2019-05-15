//------------------------------------------------------------------------------
// Look for any points that touch a surface (multi-material or void).
// For such points:
//   - Remove any non-surface multi-material overlap.
//   - If not a surface point, flag this point as touching a surface point with
//     surfacePoint=-1.
//------------------------------------------------------------------------------
#include "editMultimaterialSurfaceTopology.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"

#include <vector>

extern Timer TIME_CRKSPH_editMultimaterialSurfaceTopology;

using std::vector;

namespace Spheral {

template<typename Dimension>
void
editMultimaterialSurfaceTopology(FieldList<Dimension, int>& surfacePoint,
                                 ConnectivityMap<Dimension>& connectivityMap) {
  TIME_CRKSPH_editMultimaterialSurfaceTopology.start();

  // Declare some useful stuff and preconditions.
  const auto  numNodeLists = surfacePoint.size();
  const auto& nodeLists = connectivityMap.nodeLists();
  REQUIRE(nodeLists.size() == numNodeLists);

  // Record the points we want to remove from the neighbor set.
  FieldList<Dimension, vector<vector<int>>> neighborsToCut(FieldStorageType::CopyFields);
  for (auto nodeList: nodeLists) neighborsToCut.appendNewField("cut neighbors", *nodeList, vector<vector<int>>(numNodeLists));

  // Here we go.
  for (auto iNodeList = 0; iNodeList < numNodeLists - 1; ++iNodeList) {
    const auto n = nodeLists[iNodeList]->numInternalNodes();
    for (auto i = 0; i < n; ++i) {
      const auto& allneighbors = connectivityMap.connectivityForNode(iNodeList, i);

      // printf(" --> (%d, %d) :", iNodeList, i);

      // Assuming symmetry, we only have to visit each pair once, so we can do the whole upper triangular
      // walk here.
      for (auto jNodeList = iNodeList + 1; jNodeList < numNodeLists; ++jNodeList) {

        // If point i is a boundary for NodeList j, then it always couples to j' points.
        const bool ijboundary = (surfacePoint(iNodeList, i) && (1 << (jNodeList + 1)) > 0);
        // printf(" <%d %d> ", surfacePoint(iNodeList, i), ijboundary);
        if (not ijboundary) {

          // Since i is not a boundary for jNodeList, we need to check for any neighbors of i in
          // jNodeList that are not boundaries to iNodeList.
          const auto& neighbors = allneighbors[jNodeList];
          for (const auto j: neighbors) {
            const bool jiboundary = (surfacePoint(jNodeList, j) && (1 << (iNodeList + 1)) > 0);
            // printf(" [%d, %d, %d, %d]", jNodeList, j, surfacePoint(jNodeList, j), jiboundary);
            if (not jiboundary) {

              // Disconnect this pair.
              neighborsToCut(iNodeList, i)[jNodeList].push_back(j);
              neighborsToCut(jNodeList, j)[iNodeList].push_back(i);
              // printf(" *** Clipping (%d, %d) <--> (%d, %d)\n", iNodeList, i, jNodeList, j);

            }
          }
        }
      }
      // printf("\n");
    }
  }

  // Apply our toplogical cuts to the ConnectivityMap.
  connectivityMap.removeConnectivity(neighborsToCut);

  TIME_CRKSPH_editMultimaterialSurfaceTopology.stop();
}

}
