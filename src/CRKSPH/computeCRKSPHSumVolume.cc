//------------------------------------------------------------------------------
// Compute the CRKSPH volume summation.
//------------------------------------------------------------------------------
#include "computeCRKSPHSumVolume.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;

template<typename Dimension>
void
computeCRKSPHSumVolume(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Scalar>& vol) {

  // Pre-conditions.
  const size_t numNodeLists = vol.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;


  // For our first pass compute the effective volume per point.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(vol[nodeListi]->nodeList());

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        const int firstGhostNodej = vol[nodeListj]->nodeList().firstGhostNode();
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;

          // Check if this node pair has already been calculated.
          if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                       nodeListj, j,
                                                       firstGhostNodej)) {
            const Vector& rj = position(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Scalar etai = (Hi*rij).magnitude();
            const Scalar etaj = (Hj*rij).magnitude();
            const Scalar Wi = W.kernelValue(etai, Hdeti);
            const Scalar Wj = W.kernelValue(etaj, Hdetj);

            // Sum the pair-wise contributions.
            vol(nodeListi, i) += Wi;
            vol(nodeListj, j) += Wj;
          }
        }
      }
  
      // Add the self-contribution.
      vol(nodeListi, i) += W.kernelValue(0.0, Hdeti);
      CHECK(vol(nodeListi, i) > 0.0);
      vol(nodeListi, i) = 1.0/vol(nodeListi, i);
    }
  }

}

}
}

