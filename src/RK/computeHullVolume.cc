//---------------------------------Spheral++------------------------------------
// Compute the volume per point using an inverse convex hull.
// Optionally this volume can then be clipped to the Voronoi.
//------------------------------------------------------------------------------
#include "computeHullVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/range.hh"
#include "Geometry/PolyClipperUtilities.hh"
#include "Utilities/Timer.hh"

#include <algorithm>
#include <utility>
#include <limits>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Generic (2D, 3D) method
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeHullVolume(const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                  const ConnectivityMap<Dimension >& connectivityMap,
                  const bool clipToVoronoi,
                  FieldList<Dimension, int>& surfacePoint,
                  FieldList<Dimension, typename Dimension::Scalar>& vol,
                  FieldList<Dimension, typename Dimension::FacetedVolume>& cells) {

  TIME_FUNCTION;

  // Pre-conditions
  REQUIRE(vol.size() == position.size());

  using Vector = typename Dimension::Vector;
  using FacetedVolume = typename Dimension::FacetedVolume;
  // using PCVector = typename ClippingType<Dimension>::Vector;
  // using Plane = typename ClippingType<Dimension>::Plane;
  // using PolyVolume = typename ClippingType<Dimension>::PolyVolume;

  const auto numGens = position.numNodes();
  const auto numNodeLists = position.size();
  const auto numGensGlobal = allReduce(numGens, MPI_SUM, Communicator::communicator());
  const auto returnSurface = surfacePoint.size() == numNodeLists;
  const auto returnCells = cells.size() == numNodeLists;

  if (returnSurface) surfacePoint = 0;

  if (numGensGlobal > 0) {

    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto n = position[nodeListi]->numInternalElements();

    // Do each point independently
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hinv = Hi.Inverse();

        // Build the set of inverse positions in eta space about this point (including itself as the origin)
        vector<Vector> invPositions = {Vector::zero};
        const auto& connectivity = connectivityMap.connectivityForNode(nodeListi, i);
        CHECK(connectivity.size() == numNodeLists);
        for (auto nodeListj = 0u; nodeListj < numNodeLists; ++nodeListj) {
          for (auto j: connectivity[nodeListj]) {
            const auto etaji = (position(nodeListj, j) - ri);
            invPositions.push_back(safeInv(etaji.magnitude2()) * etaji);
          }
        }

        // Compute the inverse convex hull (in 1/eta space)
        const FacetedVolume invHull(invPositions);

        // Now we can reconstruct the inner hull in proper coordinates
        auto surface = false;
        vector<Vector> verts;
        const auto& vertsInv = invHull.vertices();
        for (const auto& vinv: vertsInv) {
          const auto vimag2 = vinv.magnitude2();
          if (vimag2 < 1.0e-30) {
            verts.push_back(Vector::zero);
            surface = true;
          } else {
            verts.push_back(vinv/vimag2);
          }
        }

        // Construct the hull in normal coordinates
        FacetedVolume hull(verts);

        // Put together our return values
        if (surface) {
          if (returnSurface) surfacePoint(nodeListi, i) = 1;
        } else {
          vol(nodeListi, i) = hull.volume();
        }
        if (returnCells) cells(nodeListi, i) = hull;
      }
    }
  }
}

}
