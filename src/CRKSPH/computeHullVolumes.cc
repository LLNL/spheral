//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on convex hulls.
//------------------------------------------------------------------------------
#include "computeHullVolumes.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

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

template<typename Dimension>
void
computeHullVolumes(const ConnectivityMap<Dimension>& connectivityMap,
                   const typename Dimension::Scalar kernelExtent,
                   const FieldList<Dimension, typename Dimension::Vector>& position,
                   const FieldList<Dimension, typename Dimension::SymTensor>& H,
                   FieldList<Dimension, typename Dimension::Scalar>& volume) {

  // Pre-conditions.
  const size_t numNodeLists = volume.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(kernelExtent > 0.0);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Collect the half-way positions of all neighbors *within i's sampling volume*.
      // We do this in eta space.
      vector<Vector> etaInv(1, Vector::zero);
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;
          const Vector& rj = position(nodeListj, j);
          const Vector rji = 0.5*(rj - ri);
          const Vector etai = Hi*rji;
          const Scalar etaiMag = etai.magnitude();
          if (etaiMag < kernelExtent) {
            const Vector etaiHat = etai.unitVector();
            etaInv.push_back(1.0/max(etaiMag, 1.0e-30) * etaiHat);
          }
        }
      }

      // Build the hull of the inverse.
      const FacetedVolume hullInv(etaInv);

      // Use the vertices selected by the inverse hull to construct the
      // volume of the node.
      vector<Vector> eta;
      const vector<Vector>& vertsInv = hullInv.vertices();
      CHECK((Dimension::nDim == 1 and vertsInv.size() == 2) or
            (Dimension::nDim == 2 and vertsInv.size() >= 3) or
            (Dimension::nDim == 3 and vertsInv.size() >= 4));
      for (typename std::vector<Vector>::const_iterator itr = vertsInv.begin();
           itr != vertsInv.end();
           ++itr) {
        if (itr->magnitude2() < 1.0e-30) {
          eta.push_back(Vector::zero);
        } else {
          eta.push_back(1.0/sqrt(itr->magnitude2()) * itr->unitVector());
        }
      }

      // And we have it.
      const FacetedVolume polyeta = FacetedVolume(eta);
      volume(nodeListi, i) = polyeta.volume()/Hdeti;
    }
  }
}

}

