//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH gradient.
//------------------------------------------------------------------------------
#include "gradientCRKSPH.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "SPH/NodeCoupling.hh"
#include "CRKSPHUtilities.hh"

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

template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientCRKSPH(const FieldList<Dimension, DataType>& fieldList,
               const FieldList<Dimension, typename Dimension::Vector>& position,
               const FieldList<Dimension, typename Dimension::Scalar>& weight,
               const FieldList<Dimension, typename Dimension::SymTensor>& H,
               const FieldList<Dimension, typename Dimension::Scalar>& A,
               const FieldList<Dimension, typename Dimension::Vector>& B,
               const FieldList<Dimension, typename Dimension::Tensor>& C,
               const FieldList<Dimension, typename Dimension::Vector>& gradA,
               const FieldList<Dimension, typename Dimension::Tensor>& gradB,
               const FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
               const ConnectivityMap<Dimension>& connectivityMap,
               const RKOrder correctionOrder,
               const TableKernel<Dimension>& W,
               const NodeCoupling& nodeCoupling) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists);
  REQUIRE(C.size() == numNodeLists or correctionOrder != RKOrder::QuadraticOrder);
  REQUIRE(gradA.size() == numNodeLists);
  REQUIRE(gradB.size() == numNodeLists or correctionOrder == RKOrder::ZerothOrder);
  REQUIRE(gradC.size() == numNodeLists or correctionOrder != RKOrder::QuadraticOrder);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Prepare the result.
  FieldList<Dimension, GradientType> result;
  result.copyFields();
  for (auto fieldItr = fieldList.begin();
       fieldItr != fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, GradientType>("grad ", (*fieldItr)->nodeList()));
  }

  // Walk the FluidNodeLists.
  auto Bi = Vector::zero, Bj = Vector::zero;
  auto Ci = Tensor::zero, Cj = Tensor::zero;
  auto gradBi = Tensor::zero, gradBj = Tensor::zero;
  auto gradCi = ThirdRankTensor::zero, gradCj = ThirdRankTensor::zero;
  for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const auto firstGhostNodei = A[nodeListi]->nodeList().firstGhostNode();

    // Iterate over the nodes in this node list.
    for (auto iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const auto i = *iItr;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto& Ai = A(nodeListi, i);
      const auto& gradAi = gradA(nodeListi, i);
      if (correctionOrder != RKOrder::ZerothOrder) {
        Bi = B(nodeListi, i);
        gradBi = gradB(nodeListi, i);
      }
      if (correctionOrder == RKOrder::QuadraticOrder) {
        Ci = C(nodeListi, i);
        gradCi = gradC(nodeListi, i);
      }
      const auto& Fi = fieldList(nodeListi, i);
      auto& gradFi = result(nodeListi, i);

      // Add our self-contribution.  A strange thing in a gradient!
      const auto W0 = W.kernelValue(0.0, Hdeti);
      gradFi += weight(nodeListi, i)*Fi*W0*(Ai*Bi + gradAi);

      // Neighbors!
      const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Walk the neighbor nodeLists.
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      
        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const auto firstGhostNodej = A[nodeListj]->nodeList().firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const auto j = *jItr;

            // The coupling between these nodes.
            const auto fij = nodeCoupling(nodeListi, i, nodeListj, j);
            
            // Only proceed if this node pair has not been calculated yet.
            if (fij > 0.0 and connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                                       nodeListj, j,
                                                                       firstGhostNodej)) {

              // Find the effective weights of i->j and j->i.
              // const Scalar wi = fij*2.0*weight(nodeListi, i)*weight(nodeListj, j)/(weight(nodeListi, i) + weight(nodeListj, j));
              // const Scalar wi = fij*0.5*(weight(nodeListi, i) + weight(nodeListj, j));
              // const Scalar wj = wi;
              const Scalar wi = fij*weight(nodeListi, i);
              const Scalar wj = fij*weight(nodeListj, j);

              // Get the state for node j.
              const auto& rj = position(nodeListj, j);
              const auto& Hj = H(nodeListj, j);
              const auto  Hdetj = Hj.Determinant();
              const auto  Aj = A(nodeListj, j);
              const auto& gradAj = gradA(nodeListj, j);
              if (correctionOrder != RKOrder::ZerothOrder) {
                Bj = B(nodeListj, j);
                gradBj = gradB(nodeListj, j);
              }
              if (correctionOrder == RKOrder::QuadraticOrder) {
                Cj = C(nodeListj, j);
                gradCj = gradC(nodeListj, j);
              }
              const auto& Fj = fieldList(nodeListj, j);
              auto& gradFj = result(nodeListj, j);

              // Node displacement.
              const auto rij = ri - rj;
              const auto etai = Hi*rij;
              const auto etaj = Hj*rij;

              // Kernel weight and gradient.
              Scalar Wi, gWi, Wj, gWj;
              Vector gradWi, gradWj;
              CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, correctionOrder,  rij,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi);
              CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, correctionOrder, -rij, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj);

              // Increment the pair-wise gradients.
              gradFi += wj*Fj*gradWj;
              gradFj += wi*Fi*gradWi;

            }
          }
        }
      }
    }
  }

  // That's it!
  return result;
}

}
