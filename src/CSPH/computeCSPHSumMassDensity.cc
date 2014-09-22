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

  // Compute the initial volume guess and prepare a FieldList to hold the zeroth correction.
  const FieldList<Dimension, Scalar> volume = mass/massDensity;
  FieldList<Dimension, Scalar> A0(FieldSpace::Copy), Veff(FieldSpace::Copy);

  // Build the correction for the initial volume estimate.
  const Scalar W0 = W.kernelValue(0.0, 1.0);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = massDensity[nodeListi]->nodeList();
    A0.appendNewField("A0", nodeList, 0.0);
    Veff.appendNewField("effective volume", nodeList, 0.0);
    const int firstGhostNodei = nodeList.firstGhostNode();
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;
      const Vector& ri = position(nodeListi, i);
      const Scalar& Vi = volume(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const vector<int>& connectivity = connectivityMap.connectivityForNode(nodeListi, i)[nodeListi];
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
        const int j = *jItr;
        if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                     nodeListi, j,
                                                     firstGhostNodei)) {
          const Vector& rj = position(nodeListi, j);
          const Scalar Vj = volume(nodeListi, j);
          const SymTensor& Hj = H(nodeListi, j);
          const Scalar Hdetj = Hj.Determinant();
          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
          const Scalar Wi = W.kernelValue(etai, Hdeti);
          const Scalar Wj = W.kernelValue(etaj, Hdetj);
          A0(nodeListi, i) += Vj*Wi;
          A0(nodeListi, j) += Vi*Wj;
        }
      }
      A0(nodeListi, i) += Vi*Hdeti*W0;
      CHECK(A0(nodeListi, i) > 0.0);
      A0(nodeListi, i) = 1.0/A0(nodeListi, i);
      // if (std::abs(A0(nodeListi, i) - 1.0) < 0.2) A0(nodeListi, i) = 1.0;
    }
  }

  // Apply boundaries to the zeroth correction.  We assume the caller has taken care of the input fields.
  for (ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(A0);
  for (ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Walk the FluidNodeLists and sum the new mass density.
  massDensity = 0.0;
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
      const Scalar Vi = volume(nodeListi, i);
      const Scalar mi = mass(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar A0i = A0(nodeListi, i);

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
          const Scalar Vj = volume(nodeListi, j);
          const Scalar mj = mass(nodeListi, j);
          const SymTensor& Hj = H(nodeListi, j);
          const Scalar Hdetj = Hj.Determinant();
          const Scalar A0j = A0(nodeListi, j);

          // Kernel weighting and gradient.
          const Vector rij = ri - rj;
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar etaj = (Hj*rij).magnitude();
          const Scalar Wi = A0i*W.kernelValue(etai, Hdeti);
          const Scalar Wj = A0j*W.kernelValue(etaj, Hdetj);

          // Sum the pair-wise contributions.
          // massDensity(nodeListi, i) += mj*Wi;
          // massDensity(nodeListi, j) += mi*Wj;
          massDensity(nodeListi, i) += Vj*mj*Wi;
          massDensity(nodeListi, j) += Vi*mi*Wj;
          Veff(nodeListi, i) += Vj*Vj*Wi;
          Veff(nodeListi, j) += Vi*Vi*Wj;
        }
      }
      
      // Finalize the density for node i.
      // massDensity(nodeListi, i) = max(rhoMin, 
      //                                 min(rhoMax,
      //                                     mi/(massDensity(nodeListi, i) + Vi*Vi*A0i*Hdeti*W0)));
      // massDensity(nodeListi, i) = max(rhoMin, 
      //                                 min(rhoMax,
      //                                     (massDensity(nodeListi, i) + mi*A0i*Hdeti*W0)));
      massDensity(nodeListi, i) = max(rhoMin, 
                                      min(rhoMax,
                                          (massDensity(nodeListi, i) + Vi*mi*A0i*Hdeti*W0) * 
                                          safeInv(Veff(nodeListi, i) + Vi*Vi*A0i*Hdeti*W0)));
      CHECK(massDensity(nodeListi, i) > 0.0);
    }
  }
}

}
}

