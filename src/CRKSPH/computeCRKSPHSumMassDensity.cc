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
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  // // Compute the CRK corrections.
  // FieldList<Dimension, Scalar>           m0(FieldSpace::Copy);
  // FieldList<Dimension, Vector>           m1(FieldSpace::Copy);
  // FieldList<Dimension, SymTensor>        m2(FieldSpace::Copy);
  // FieldList<Dimension, ThirdRankTensor>  m3(FieldSpace::Copy);
  // FieldList<Dimension, FourthRankTensor> m4(FieldSpace::Copy);
  // FieldList<Dimension, Vector>           gradm0(FieldSpace::Copy);
  // FieldList<Dimension, Tensor>           gradm1(FieldSpace::Copy);
  // FieldList<Dimension, ThirdRankTensor>  gradm2(FieldSpace::Copy);
  // FieldList<Dimension, FourthRankTensor> gradm3(FieldSpace::Copy);
  // FieldList<Dimension, FifthRankTensor>  gradm4(FieldSpace::Copy);
  // FieldList<Dimension, Scalar>           A(FieldSpace::Copy);
  // FieldList<Dimension, Vector>           B(FieldSpace::Copy);
  // FieldList<Dimension, Tensor>           C(FieldSpace::Copy);
  // FieldList<Dimension, Vector>           gradA(FieldSpace::Copy);
  // FieldList<Dimension, Tensor>           gradB(FieldSpace::Copy);
  // FieldList<Dimension, ThirdRankTensor>  gradC(FieldSpace::Copy);
  // FieldList<Dimension, Scalar>           one(FieldSpace::Copy);
  // for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //   m0.appendNewField(HydroFieldNames::m0_CRKSPH, position[nodeListi]->nodeList(), 0.0);
  //   m1.appendNewField(HydroFieldNames::m1_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   m2.appendNewField(HydroFieldNames::m2_CRKSPH, position[nodeListi]->nodeList(), SymTensor::zero);
  //   m3.appendNewField(HydroFieldNames::m3_CRKSPH, position[nodeListi]->nodeList(), ThirdRankTensor::zero);
  //   m4.appendNewField(HydroFieldNames::m4_CRKSPH, position[nodeListi]->nodeList(), FourthRankTensor::zero);
  //   gradm0.appendNewField(HydroFieldNames::m0_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   gradm1.appendNewField(HydroFieldNames::m1_CRKSPH, position[nodeListi]->nodeList(), Tensor::zero);
  //   gradm2.appendNewField(HydroFieldNames::m2_CRKSPH, position[nodeListi]->nodeList(), ThirdRankTensor::zero);
  //   gradm3.appendNewField(HydroFieldNames::m3_CRKSPH, position[nodeListi]->nodeList(), FourthRankTensor::zero);
  //   gradm4.appendNewField(HydroFieldNames::m4_CRKSPH, position[nodeListi]->nodeList(), FifthRankTensor::zero);
  //   A.appendNewField(HydroFieldNames::A_CRKSPH, position[nodeListi]->nodeList(), 0.0);
  //   B.appendNewField(HydroFieldNames::B_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   C.appendNewField(HydroFieldNames::C_CRKSPH, position[nodeListi]->nodeList(), Tensor::zero);
  //   gradA.appendNewField(HydroFieldNames::gradA_CRKSPH, position[nodeListi]->nodeList(), Vector::zero);
  //   gradB.appendNewField(HydroFieldNames::gradB_CRKSPH, position[nodeListi]->nodeList(), Tensor::zero);
  //   gradC.appendNewField(HydroFieldNames::gradC_CRKSPH, position[nodeListi]->nodeList(), ThirdRankTensor::zero);
  //   one.appendNewField("one", position[nodeListi]->nodeList(), 1.0);
  // }
  // const FieldList<Dimension, Scalar> vol0 = mass/massDensity;
  // PerNodeListNodeCoupling coupling;
  // computeCRKSPHMoments(connectivityMap, W, vol0, position, H, QuadraticOrder, coupling,
  //                      m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
  // computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, QuadraticOrder,
  //                          A, B, C, gradA, gradB, gradC);

  // massDensity = interpolateCRKSPH(mass, position, one, H, A, B, C, connectivityMap, QuadraticOrder, W, coupling);

  // const FieldList<Dimension, Scalar> mavg = interpolateCRKSPH(mass, position, vol0, H, A, B, C, connectivityMap, ZerothOrder, W, coupling);
  // const FieldList<Dimension, Scalar> volavg = interpolateCRKSPH(vol0, position, mass, H, A, B, C, connectivityMap, ZerothOrder, W, coupling);

  // for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //   const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
  //   const Scalar rhoMin = nodeList.rhoMin();
  //   const Scalar rhoMax = nodeList.rhoMax();
  //   const size_t n = nodeList.numInternalNodes();
  //   for (size_t i = 0; i != n; ++i) {
  //     CHECK(volavg(nodeListi,i) > 0.0);
  //     massDensity(nodeListi, i) = max(rhoMin, min(rhoMax, mavg(nodeListi, i))); // /volavg(nodeListi, i)));
  //   }
  // }

  FieldList<Dimension, Scalar> vol0(FieldSpace::Copy);
  FieldList<Dimension, Scalar> vol1(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    vol0.appendNewField("volume", position[nodeListi]->nodeList(), 0.0);
    vol1.appendNewField("averaged volume", position[nodeListi]->nodeList(), 0.0);
  }

  // For our first pass compute the effective volume per point.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());

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
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Scalar etai = (Hi*rij).magnitude();
            const Scalar etaj = (Hj*rij).magnitude();
            const Scalar Wi = W.kernelValue(etai, Hdeti);
            const Scalar Wj = W.kernelValue(etaj, Hdetj);

            // Sum the pair-wise contributions.
            vol0(nodeListi, i) += Wi;
            vol0(nodeListj, j) += Wj;
          }
        }
      }
  
      // Add the self-contribution.
      vol0(nodeListi, i) += W.kernelValue(0.0, Hdeti);
      CHECK(vol0(nodeListi, i) > 0.0);
      vol0(nodeListi, i) = 1.0/vol0(nodeListi, i);
    }
  }

  // Apply boundary conditions to the per node volume estimate.
  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator BoundaryIterator;
  for (BoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(vol0);
  for (BoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Now use this volume as weighting to compute the effective volume and mass per point.
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
      const Scalar Vi = vol0(nodeListi, i);
      const Scalar rhoi = massDensity(nodeListi, i);
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
            const Scalar Vj = vol0(nodeListj, j);
            const Scalar rhoj = massDensity(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            // Kernel weighting and gradient.
            const Vector rij = ri - rj;
            const Scalar etai = (Hi*rij).magnitude();
            const Scalar etaj = (Hj*rij).magnitude();
            const Scalar Wi = W.kernelValue(etai, Hdeti);
            const Scalar Wj = W.kernelValue(etaj, Hdetj);

            // Sum the pair-wise contributions.
            massDensity(nodeListi, i) += (nodeListi == nodeListj ? mj : mi) * Vj*Wi;
            massDensity(nodeListj, j) += (nodeListi == nodeListj ? mi : mj) * Vi*Wj;
            vol1(nodeListi, i) += Vj*Vj*Wi;
            vol1(nodeListj, j) += Vi*Vi*Wj;
          }
        }
      }
  
      // Finalize the density for node i.
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + mi*Vi*W.kernelValue(0.0, Hdeti))/
                                          (vol1(nodeListi, i) + Vi*Vi*W.kernelValue(0.0, Hdeti))));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
}

