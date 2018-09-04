//------------------------------------------------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#include "computeSolidCRKSPHSumMassDensity.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "SPH/NodeCoupling.hh"

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
computeSolidCRKSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                                 const TableKernel<Dimension>& W,
                                 const FieldList<Dimension, typename Dimension::Vector>& position,
                                 const FieldList<Dimension, typename Dimension::Scalar>& mass,
                                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                                 const FieldList<Dimension, typename Dimension::Scalar>& massDensity0,
                                 const NodeCoupling& nodeCoupling,
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
  typedef typename std::vector<Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  const Scalar W0 = W.kernelValue(0.0, 1.0);

  // Prepare to sum the correction.
  FieldList<Dimension, Scalar> m0(FieldStorageType::CopyFields);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    m0.appendNewField("zeroth correction", position[nodeListi]->nodeList(), 0.0);
  }

  // Walk the FluidNodeLists and sum the new mass density.
  massDensity = 0.0;
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
      const Scalar rho0i = massDensity0(nodeListi, i);
      const Scalar wi = mi/rho0i;
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        const int firstGhostNodej = massDensity[nodeListj]->nodeList().firstGhostNode();
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;

          // Check the coupling of these points.
          const Scalar fij = nodeCoupling(nodeListi, i, nodeListj, j);

          // Check if this node pair has already been calculated.
          if (fij > 0.0 and connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                                     nodeListj, j,
                                                                     firstGhostNodej)) {
            const Vector& rj = position(nodeListj, j);
            const Scalar mj = mass(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();
            const Scalar rho0j = massDensity0(nodeListj, j);
            const Scalar wj = mj/rho0j;

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Scalar etai = (Hi*rij).magnitude();
            const Scalar etaj = (Hj*rij).magnitude();
            const Scalar Wi = W.kernelValue(etai, Hdeti);
            const Scalar Wj = W.kernelValue(etaj, Hdetj);

            // Sum the pair-wise contributions.
            massDensity(nodeListi, i) += fij*wj*Wj*mi;
            massDensity(nodeListj, j) += fij*wi*Wi*mj;
            m0(nodeListi, i) += fij*wj*Wj*wi;
            m0(nodeListj, j) += fij*wi*Wi*wj;
          }
        }
      }
      
      // Finalize the density for node i.
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + wi*Hdeti*W0*mi)/
                                          (m0(nodeListi, i) + wi*Hdeti*W0*wi)));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}

