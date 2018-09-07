//------------------------------------------------------------------------------
// Compute the TaylorSPH corrections.
//------------------------------------------------------------------------------
#include "computeTaylorSPHCorrections.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerProduct.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

template<typename Dimension>
void
computeTaylorSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Scalar>& weight,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Tensor>& D) {

  // Pre-conditions.
  const size_t numNodeLists = D.size();
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Zero out the result.
  D = Tensor::zero;

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Scalar wi = weight(nodeListi, i);
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Neighbors!
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const NodeList<Dimension>& nodeList = D[nodeListj]->nodeList();
        const int firstGhostNodej = nodeList.firstGhostNode();
        const vector<int>& connectivity = fullConnectivity[nodeListj];

        // Iterate over the neighbors for in this NodeList.
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;

          // Check if this node pair has already been calculated.
          if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                       nodeListi, j,
                                                       firstGhostNodej)) {

            // State of node j.
            const Scalar wj = weight(nodeListj, j);
            const Vector& rj = position(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;
            const std::pair<double, double> WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
            const Scalar& Wi = WWi.first;
            const Vector gradWi = (Hi*etai.unitVector())*WWi.second;
            const std::pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
            const Scalar& Wj = WWj.first;
            const Vector gradWj = (Hj*etaj.unitVector())*WWj.second;

            D(nodeListi, i) -= wj*rij.dyad(gradWj);
            D(nodeListj, j) -= wi*rij.dyad(gradWi);
          }
        }
      }
      D(nodeListi, i) = D(nodeListi, i).Inverse();
    }
  }
}

}

