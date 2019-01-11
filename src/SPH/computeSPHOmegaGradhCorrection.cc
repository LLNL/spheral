//------------------------------------------------------------------------------
// Compute the SPH grad h correction due to Springel et al.
//------------------------------------------------------------------------------
#include "computeSPHOmegaGradhCorrection.hh"
#include "Field/Field.hh"
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
computeSPHOmegaGradhCorrection(const ConnectivityMap<Dimension>& connectivityMap,
                               const TableKernel<Dimension>& W,
                               const FieldList<Dimension, typename Dimension::Vector>& position,
                               const FieldList<Dimension, typename Dimension::SymTensor>& H,
                               FieldList<Dimension, typename Dimension::Scalar>& omegaGradh) {

  // Pre-conditions.
  const size_t numNodeLists = omegaGradh.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);
  Scalar Wi, gWi;

  // Zero out the result.
  omegaGradh = 0.0;

  // Prepare a FieldList to hold the sum gradient.
  FieldList<Dimension, Scalar> gradsum(FieldStorageType::CopyFields);
  for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const auto& nodeList = omegaGradh[nodeListi]->nodeList();
    gradsum.appendNewField("sum of the gradient", nodeList, 0.0);
  }

  // Walk the FluidNodeLists.
  for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const auto& nodeList = omegaGradh[nodeListi]->nodeList();
    const auto ni = connectivityMap.numNodes(nodeListi);

    // Iterate over the nodes in this node list.
#pragma omp parallel for private(Wi, gWi)
    for (auto k = 0; k < ni; ++k) {
      const auto i = connectivityMap.ithNode(nodeListi, k);

      // If we're isolated we have to punt and just set this correction to unity.
      if (connectivityMap.numNeighborsForNode(nodeListi, i) == 0) {
        omegaGradh(nodeListi, i) = 1.0;

      } else {

        // Get the state for node i.
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();

        // Self-contribution.
        omegaGradh(nodeListi, i) += Hdeti*W0;

        // Neighbors!
        const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        CHECK(fullConnectivity.size() == numNodeLists);

        // Iterate over the neighbor NodeLists.
        for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          const auto firstGhostNodej = omegaGradh[nodeListj]->nodeList().firstGhostNode();

          // Iterate over the neighbors for in this NodeList.
          const auto& connectivity = fullConnectivity[nodeListj];
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const auto j = *jItr;
            const auto& rj = position(nodeListj, j);

            // Kernel weighting and gradient.
            const auto rij = ri - rj;
            const auto etai = (Hi*rij).magnitude();
            std::tie(Wi, gWi) = W.kernelAndGradValue(etai, Hdeti);

            // Sum the pair-wise contributions.
            omegaGradh(nodeListi, i) += Wi;
            gradsum(nodeListi, i) += etai*gWi;
          }
        }

        // Finish the grad h correction.
        CHECK(omegaGradh(nodeListi, i) > 0.0);
        omegaGradh(nodeListi, i) = std::max(1.0e-30, -gradsum(nodeListi, i)/(Dimension::nDim * omegaGradh(nodeListi, i)));
      }

      // Post-conditions.
      ENSURE2(i >= nodeList.firstGhostNode() or omegaGradh(nodeListi, i) >= 0.0, 
              nodeListi << " " << i << " " << nodeList.firstGhostNode() << " " << omegaGradh(nodeListi, i));
    }
  }
}

}
