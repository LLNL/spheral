//------------------------------------------------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#include "computeCRKSPHSumMassDensity.hh"
#include "interpolateCRKSPH.hh"
#include "computeCRKSPHMoments.hh"
#include "computeCRKSPHCorrections.hh"
#include "CRKSPHUtilities.hh"
#include "CRKSPHCorrectionParams.hh"
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
                            const FieldList<Dimension, typename Dimension::Scalar>& vol,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            const FieldList<Dimension, typename Dimension::Scalar>& A,
                            const FieldList<Dimension, typename Dimension::Vector>& B,
                            const FieldList<Dimension, typename Dimension::Tensor>& C,
                            const CRKOrder correctionOrder,
                            FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(vol.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(A.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists or correctionOrder == ZerothOrder);
  REQUIRE(C.size() == numNodeLists or correctionOrder != QuadraticOrder);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  FieldList<Dimension, Scalar> vol1(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    vol1.appendNewField("averaged volume", position[nodeListi]->nodeList(), 0.0);
  }

  // Some scratch variables.
  Scalar Wi, Wj;
  Vector rij, etai, etaj;
  Vector Bi = Vector::zero, Bj = Vector::zero;
  Tensor Ci = Tensor::zero, Cj = Tensor::zero;

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
      const Scalar Vi = vol(nodeListi, i);
      const Scalar rhoi = massDensity(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar Ai = A(nodeListi, i);
      if (correctionOrder != ZerothOrder) {
        Bi = B(nodeListi, i);
      }
      if (correctionOrder == QuadraticOrder) {
        Ci = C(nodeListi, i);
      }
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
            const Scalar Vj = vol(nodeListj, j);
            const Scalar rhoj = massDensity(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();
            const Scalar Aj = A(nodeListj, j);
            if (correctionOrder != ZerothOrder) {
              Bj = B(nodeListj, j);
            }
            if (correctionOrder == QuadraticOrder) {
              Cj = C(nodeListj, j);
            }

            // Kernel weighting and gradient.
            rij = ri - rj;
            etai = Hi*rij;
            etaj = Hj*rij;
            // Wj = CRKSPHKernel(W, correctionOrder,  rij,  etai, Hdeti,  etaj, Hdetj, Ai, Bi, Ci);
            // Wi = CRKSPHKernel(W, correctionOrder, -rij, -etaj, Hdetj, -etai, Hdeti, Aj, Bj, Cj);
            Wj = W.kernelValue(etaj.magnitude(), Hdetj);
            Wi = W.kernelValue(etai.magnitude(), Hdeti);

            // Sum the pair-wise contributions.
            massDensity(nodeListi, i) += (nodeListi == nodeListj ? mj : mi) * Vj*Wj;
            massDensity(nodeListj, j) += (nodeListi == nodeListj ? mi : mj) * Vi*Wi;
            vol1(nodeListi, i) += Vj*Vj*Wj;
            vol1(nodeListj, j) += Vi*Vi*Wi;
          }
        }
      }
  
      // Finalize the density for node i.
      // Wj = CRKSPHKernel(W, correctionOrder, Vector::zero, Vector::zero, Hdeti, Vector::zero, Hdeti, Ai, Bi, Ci);
      Wj = W.kernelValue(0.0, Hdeti);
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + mi*Vi*Wj)/
                                          (vol1(nodeListi, i) + Vi*Vi*Wj)));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
}

