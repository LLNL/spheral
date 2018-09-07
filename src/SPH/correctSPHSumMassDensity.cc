//------------------------------------------------------------------------------
// Compute the SPH mass density summation.
//------------------------------------------------------------------------------
#include "computeSPHSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
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
correctSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const bool sumOverAllNodeLists,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::Scalar>& mass,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Make a single corrective pass.
  FieldList<Dimension, Scalar> sumUnity(FieldStorageType::CopyFields);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    sumUnity.appendNewField("SPH sum unity check", massDensity[nodeListi]->nodeList(), 0.0);
  }
  
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = massDensity[nodeListi]->nodeList();
    const int firstGhostNodei = nodeList.firstGhostNode();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar mi = mass(nodeListi, i);
      const Scalar rhoi = massDensity(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      CHECK(rhoi > 0.0);

      // Self-contribution.
      const Scalar W0 = W.kernelValue(0.0, Hdeti);
      sumUnity(nodeListi, i) += mi/rhoi*W0;

      // Get the neighbors for this node.
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        if (sumOverAllNodeLists or (nodeListi == nodeListj)) {
          const int firstGhostNodej = massDensity[nodeListj]->nodeList().firstGhostNode();
          const vector<int>& connectivity = fullConnectivity[nodeListj];
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Check if this node pair has already been calculated.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              const Vector& rj = position(nodeListj, j);
              const Scalar mj = mass(nodeListj, j);
              const Scalar rhoj = massDensity(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              CHECK(rhoj > 0.0);

              // Kernel weighting and gradient.
              const Vector rij = ri - rj;
              const Scalar etai = (Hi*rij).magnitude();
              const Scalar etaj = (Hj*rij).magnitude();
              const Scalar Wi = W.kernelValue(etai, Hdeti);
              const Scalar Wj = W.kernelValue(etaj, Hdetj);

              // Sum the pair-wise contributions.
              sumUnity(nodeListi, i) += mj/rhoj*Wj;
              sumUnity(nodeListj, j) += mi/rhoi*Wi;
            }
          }
        }
      }
      CHECK(sumUnity(nodeListi, i) > 0.0);
    }
  }

  // Apply the correction.
  massDensity /= sumUnity;
}

}
