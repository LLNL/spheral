//------------------------------------------------------------------------------
// Compute the SPH mass density summation.
//------------------------------------------------------------------------------

#include "computeSPHSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

template<typename Dimension>
void
computeSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                         const TableKernel<Dimension>& W,
                         const bool sumOverAllNodeLists,
                         const FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldList<Dimension, typename Dimension::Scalar>& mass,
                         const FieldList<Dimension, typename Dimension::SymTensor>& H,
                         FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const auto numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);

  // Zero out the result.
  massDensity = 0.0;

  // Walk the FluidNodeLists.
  for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const auto& nodeList = massDensity[nodeListi]->nodeList();
    const auto  ni = connectivityMap.numNodes(nodeListi);

    // Iterate over the nodes in this node list.
#pragma omp parallel for
    for (auto k = 0; k < ni; ++k) {
      const auto i = connectivityMap.ithNode(nodeListi, k);

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // Self-contribution.
      massDensity(nodeListi, i) += mi*Hdeti*W0;

      // Get the neighbors for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        if (sumOverAllNodeLists or (nodeListi == nodeListj)) {
          const auto& connectivity = fullConnectivity[nodeListj];
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const auto j = *jItr;

            const auto& rj = position(nodeListj, j);
            const auto  mj = mass(nodeListj, j);
            const auto& Hj = H(nodeListj, j);
            const auto  Hdetj = Hj.Determinant();

            // Kernel weighting and gradient.
            const auto rij = ri - rj;
            const auto etaj = (Hj*rij).magnitude();
            const auto Wj = W.kernelValue(etaj, Hdetj);

            // Sum the pair-wise contributions.
            massDensity(nodeListi, i) += (nodeListi == nodeListj ? mj : mi)*Wj;
          }
        }
      }
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
