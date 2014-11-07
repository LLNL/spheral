//------------------------------------------------------------------------------
// Compute the CSPH mass density summation.
//------------------------------------------------------------------------------

#include "computeCSPHSumMassDensity.hh"
#include "computeCSPHCorrections.hh"
#include "CSPHUtilities.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {
namespace CSPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using BoundarySpace::Boundary;

template<typename Dimension>
void
computeCSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                          const TableKernel<Dimension>& W,
                          const FieldList<Dimension, typename Dimension::Vector>& position,
                          const FieldList<Dimension, typename Dimension::Scalar>& mass,
                          const FieldList<Dimension, typename Dimension::SymTensor>& H,
                          const typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator& boundaryBegin,
                          const typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator& boundaryEnd,
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
  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  massDensity = 0.0;
  const Scalar W0 = W.kernelValue(0.0, 1.0);

  // Walk the FluidNodeLists and sum the new mass density.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
    const Scalar rhoMin = nodeList.rhoMin();
    const Scalar rhoMax = nodeList.rhoMax();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar mi = mass(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        const int firstGhostNodej = massDensity[nodeListj]->nodeList().firstGhostNode();
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
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Scalar etai = (Hi*rij).magnitude();
            const Scalar etaj = (Hj*rij).magnitude();
            const Scalar Wi = W.kernelValue(etai, Hdeti);
            const Scalar Wj = W.kernelValue(etaj, Hdetj);

            // Sum the pair-wise contributions.
            massDensity(nodeListi, i) += Wi;
            massDensity(nodeListj, j) += Wj;
            // massDensity(nodeListi, i) += mj*(Wi + Wj);
            // massDensity(nodeListj, j) += mi*(Wi + Wj);
            // Veff(nodeListi, i) += Vj*(Wi + Wj);
            // Veff(nodeListj, j) += Vi*(Wi + Wj);
          }
        }
      }
      
      // Finalize the density for node i.
      // massDensity(nodeListi, i) = max(rhoMin, 
      //                                 min(rhoMax,
      //                                     mi/(massDensity(nodeListi, i) + Vi*Vi*A0i*Hdeti*W0)));
      // massDensity(nodeListi, i) = max(rhoMin, 
      //                                 min(rhoMax,
      //                                     A0i*(massDensity(nodeListi, i) + mi*2.0*Hdeti*W0)));
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          mi*(massDensity(nodeListi, i) + Hdeti*W0)));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
}

