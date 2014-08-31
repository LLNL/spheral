//------------------------------------------------------------------------------
// Compute the CSPH mass density summation.
//------------------------------------------------------------------------------

#include "computeCSPHSumMassDensity.hh"
#include "computeCSPHCorrections.hh"
#include "CSPHUtilities.hh"
#include "interpolateCSPH.hh"
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
                          const FieldList<Dimension, typename Dimension::Scalar>& volume0,
                          const FieldList<Dimension, typename Dimension::SymTensor>& H,
                          const typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator& boundaryBegin,
                          const typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator& boundaryEnd,
                          FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const size_t numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(volume0.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename std::vector<BoundarySpace::Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  // Zero out the result, and prepare a FieldList to hold the effective volume.
  massDensity = 0.0;
  FieldList<Dimension, Scalar> Veff(FieldSpace::Copy);
  FieldList<Dimension, Scalar> m0(FieldSpace::Copy);
  FieldList<Dimension, Vector> m1(FieldSpace::Copy);
  FieldList<Dimension, SymTensor> m2(FieldSpace::Copy);
  FieldList<Dimension, Scalar> A0(FieldSpace::Copy);
  FieldList<Dimension, Scalar> A(FieldSpace::Copy);
  FieldList<Dimension, Vector> B(FieldSpace::Copy);
  FieldList<Dimension, Vector> C(FieldSpace::Copy);
  FieldList<Dimension, Tensor> D(FieldSpace::Copy);
  FieldList<Dimension, Vector> gradA(FieldSpace::Copy);
  FieldList<Dimension, Tensor> gradB(FieldSpace::Copy);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = massDensity[nodeListi]->nodeList();
    Veff.appendNewField("effective volume", nodeList, 0.0);
    m0.appendNewField(HydroFieldNames::m0_CSPH, nodeList, 0.0);
    m1.appendNewField(HydroFieldNames::m1_CSPH, nodeList, Vector::zero);
    m2.appendNewField(HydroFieldNames::m2_CSPH, nodeList, SymTensor::zero);
    A0.appendNewField(HydroFieldNames::A0_CSPH, nodeList, 0.0);
    A.appendNewField(HydroFieldNames::A_CSPH, nodeList, 0.0);
    B.appendNewField(HydroFieldNames::B_CSPH, nodeList, Vector::zero);
    C.appendNewField(HydroFieldNames::C_CSPH, nodeList, Vector::zero);
    D.appendNewField(HydroFieldNames::D_CSPH, nodeList, Tensor::zero);
    gradA.appendNewField(HydroFieldNames::gradA_CSPH, nodeList, Vector::zero);
    gradB.appendNewField(HydroFieldNames::gradB_CSPH, nodeList, Tensor::zero);
  }

  // First up we need to compute the CSPH corrections.  We're actually only going to use the constant 
  // correction A0.  We also force the CSPH corrections to *not* couple across NodeLists -- each NodeList
  // acts as though it is isolated.
  computeCSPHCorrections(connectivityMap, W, volume0, position, H, false,
                         m0, m1, m2, A0, A, B, C, D, gradA, gradB);

  // Apply boundaries to the zeroth correction.  We assume the caller has taken care of the input fields.
  for (ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(A0);
  for (ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // We take the smoothed volume information.
  // FieldList<Dimension, Vector> Bzero(FieldSpace::Copy);
  // for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //   const NodeList<Dimension>& nodeList = massDensity[nodeListi]->nodeList();
  //   Bzero.appendNewField("B0", nodeList, Vector::zero);
  // }
  FieldList<Dimension, Scalar> volume1 = interpolateCSPH(volume0, position, volume0, H, A, B, connectivityMap, W);

  for (ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(volume1);
  for (ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = dynamic_cast<const FluidNodeList<Dimension>&>(massDensity[nodeListi]->nodeList());
    const int firstGhostNodei = nodeList.firstGhostNode();
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
      const Scalar& V0i = volume0(nodeListi, i);
      const Scalar& V1i = volume1(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar A0i = A0(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Get the neighbors for this node (in this NodeList).  We use the approximation here
      // that nodes from other NodeLists do not contribute to the density of this one.
      const vector<int>& connectivity = connectivityMap.connectivityForNode(nodeListi, i)[nodeListi];
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
        const int j = *jItr;

        // Check if this node pair has already been calculated.
        if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                     nodeListi, j,
                                                     firstGhostNodei)) {
          const Vector& rj = position(nodeListi, j);
          const Scalar mj = mass(nodeListi, j);
          const Scalar V0j = volume0(nodeListi, j);
          const Scalar V1j = volume1(nodeListi, j);
          const SymTensor& Hj = H(nodeListi, j);
          const Scalar A0j = A0(nodeListi, j);
          const Scalar Hdetj = Hj.Determinant();

          // Kernel weighting and gradient.
          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
          const Scalar Wi = CSPHKernel(W, rij, etai, Hdeti, A0j, Vector::zero);
          const Scalar Wj = CSPHKernel(W, rij, etaj, Hdetj, A0i, Vector::zero);

          // Sum the pair-wise contributions.
          Veff(nodeListi, i) += V0j*V1j*Wj;
          massDensity(nodeListi, i) += V0j*mj*Wj;

          Veff(nodeListi, j) += V0i*V1i*Wi;
          massDensity(nodeListi, j) += V0i*mi*Wi;
        }
      }
      
      // Finalize the density for node i.
      const Scalar W0 = CSPHKernel(W, Vector::zero, 0.0, Hdeti, A0i, Vector::zero);
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + V0i*mi*W0) * 
                                          safeInv(A0i*(Veff(nodeListi, i) + V0i*V1i*W0))));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
}

