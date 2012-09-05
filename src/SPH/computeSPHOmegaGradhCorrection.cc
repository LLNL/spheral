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
namespace SPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;

template<typename Dimension>
void
computeSPHOmegaGradhCorrection(const ConnectivityMap<Dimension>& connectivityMap,
                               const TableKernel<Dimension>& W,
                               const FieldList<Dimension, typename Dimension::Vector>& position,
                               const FieldList<Dimension, typename Dimension::SymTensor>& H,
                               FieldList<Dimension, typename Dimension::Scalar>& omegaGradh) {

  // TAU timers.
  TAU_PROFILE("computeSPHOmegaGradhCorrection", "", TAU_USER);

  // Pre-conditions.
  const size_t numNodeLists = omegaGradh.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Zero out the result.
  omegaGradh = 0.0;

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = omegaGradh[nodeListi]->nodeList();

    // Stuff we're going to accumulate.
    Field<Dimension, Scalar> gradsum("sum of the gradient", nodeList);

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Self-contribution.
      const Scalar W0 = W.kernelValue(0.0, Hdeti);
      omegaGradh(nodeListi, i) += W0;

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Iterate over the neighbor NodeLists.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Iterate over the neighbors for in this NodeList.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;

          const Vector& rj = position(nodeListj, j);
          const SymTensor& Hj = H(nodeListj, j);
          const Scalar Hdetj = Hj.Determinant();

          // Kernel weighting and gradient.
          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
          const std::pair<double, double> WWi = W.kernelAndGradValue(etai, Hdeti);
          const Scalar& Wi = WWi.first;
          const Scalar& gWi = WWi.second;
          const std::pair<double, double> WWj = W.kernelAndGradValue(etaj, Hdetj);
          const Scalar& Wj = WWj.first;
          const Scalar& gWj = WWj.second;

          // Sum the pair-wise contributions.
          omegaGradh(nodeListi, i) += Wi;
          gradsum(i) += etai*gWi;
        }
      }

      // Finish the grad h correction.
      CHECK(omegaGradh(nodeListi, i) > 0.0);
      omegaGradh(nodeListi, i) = std::max(1.0e-50, -gradsum(i)/(Dimension::nDim * omegaGradh(nodeListi, i)));

      // Post-conditions.
      ENSURE(omegaGradh(nodeListi, i) > 0.0);
    }
  }
}

}
}

