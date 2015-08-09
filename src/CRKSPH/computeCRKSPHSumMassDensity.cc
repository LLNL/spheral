//------------------------------------------------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#include "computeCRKSPHSumMassDensity.hh"
#include "interpolateCRKSPH.hh"
#include "CRKSPHUtilities.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
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
computeCRKSPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                          const TableKernel<Dimension>& W,
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

  // // Compute an effective mass per point.
  // FieldList<Dimension, Scalar> m0(FieldSpace::Copy);
  // FieldList<Dimension, Vector> m1(FieldSpace::Copy);
  // FieldList<Dimension, SymTensor> m2(FieldSpace::Copy);
  // FieldList<Dimension, Scalar> A0(FieldSpace::Copy);
  // FieldList<Dimension, Scalar> A(FieldSpace::Copy);
  // FieldList<Dimension, Vector> B(FieldSpace::Copy);
  // FieldList<Dimension, Vector> C(FieldSpace::Copy);
  // FieldList<Dimension, Tensor> D(FieldSpace::Copy);
  // FieldList<Dimension, Vector> gradA0(FieldSpace::Copy);
  // FieldList<Dimension, Vector> gradA(FieldSpace::Copy);
  // FieldList<Dimension, Tensor> gradB(FieldSpace::Copy);
  // for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //   m0.appendNewField(HydroFieldNames::m0_CRKSPH, position[nodeListi]->nodeList(), 0.0);
  //   m1.appendNewField(HydroFieldNames::m1_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   m2.appendNewField(HydroFieldNames::m2_CRKSPH, position[nodeListi]->nodeList(), SymTensor::zero);
  //   A0.appendNewField(HydroFieldNames::A0_CRKSPH, position[nodeListi]->nodeList(), 0.0);
  //   A.appendNewField(HydroFieldNames::A_CRKSPH, position[nodeListi]->nodeList(), 0.0);
  //   B.appendNewField(HydroFieldNames::B_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   C.appendNewField(HydroFieldNames::C_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   D.appendNewField(HydroFieldNames::D_CRKSPH, position[nodeListi]->nodeList(), Tensor::zero);
  //   gradA0.appendNewField(HydroFieldNames::gradA0_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   gradA.appendNewField(HydroFieldNames::gradA_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   gradB.appendNewField(HydroFieldNames::gradB_CRKSPH, position[nodeListi]->nodeList(), Tensor::zero);
  // }
  // computeCRKSPHCorrections(connectivityMap, W, mass, position, H, false, m0, m1, m2,
  //                        A0, A, B, C, D, gradA0, gradA, gradB);
  // const FieldList<Dimension, Scalar> mavg = interpolateCRKSPH(mass, position, mass, H, false, A, B, connectivityMap, W);

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
            massDensity(nodeListi, i) += (nodeListi == nodeListj ? mj : mi) * Wi;
            massDensity(nodeListj, j) += (nodeListi == nodeListj ? mi : mj) * Wj;
          }
        }
      }
      
      // Finalize the density for node i.
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          massDensity(nodeListi, i) + mi*W.kernelValue(0.0, Hdeti)));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
}

