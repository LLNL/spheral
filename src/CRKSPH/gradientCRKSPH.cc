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

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;

template<typename Dimension, typename DataType>
FieldSpace::FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientCRKSPH(const FieldSpace::FieldList<Dimension, DataType>& fieldList,
               const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
               const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
               const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
               const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
               const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
               const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& C,
               const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
               const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB,
               const FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
               const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
               const CRKOrder correctionOrder,
               const KernelSpace::TableKernel<Dimension>& W,
               const NodeCoupling& nodeCoupling) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists);
  REQUIRE(C.size() == numNodeLists or correctionOrder != CRKOrder::QuadraticOrder);
  REQUIRE(gradA.size() == numNodeLists);
  REQUIRE(gradB.size() == numNodeLists or correctionOrder == CRKOrder::ZerothOrder);
  REQUIRE(gradC.size() == numNodeLists or correctionOrder != CRKOrder::QuadraticOrder);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Prepare the result.
  FieldList<Dimension, GradientType> result;
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr != fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, GradientType>("grad", (*fieldItr)->nodeList()));
  }

  // Walk the FluidNodeLists.
  Vector Bi = Vector::zero, Bj = Vector::zero;
  Tensor Ci = Tensor::zero, Cj = Tensor::zero;
  Tensor gradBi = Tensor::zero, gradBj = Tensor::zero;
  ThirdRankTensor gradCi = ThirdRankTensor::zero, gradCj = ThirdRankTensor::zero;
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const int firstGhostNodei = A[nodeListi]->nodeList().firstGhostNode();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar& Ai = A(nodeListi, i);
      const Vector& gradAi = gradA(nodeListi, i);
      if (correctionOrder != CRKOrder::ZerothOrder) {
        Bi = B(nodeListi, i);
        gradBi = gradB(nodeListi, i);
      }
      if (correctionOrder == CRKOrder::QuadraticOrder) {
        Ci = C(nodeListi, i);
        gradCi = gradC(nodeListi, i);
      }
      const DataType& Fi = fieldList(nodeListi, i);
      GradientType& gradFi = result(nodeListi, i);

      // Add our self-contribution.  A strange thing in a gradient!
      const Scalar W0 = W.kernelValue(0.0, Hdeti);
      gradFi += weight(nodeListi, i)*Fi*W0*(Ai*Bi + gradAi);

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Walk the neighbor nodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      
        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
	  const int firstGhostNodej = A[nodeListj]->nodeList().firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // The coupling between these nodes.
            const double fij = nodeCoupling(nodeListi, i, nodeListj, j);
            
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
              const Vector& rj = position(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar& Aj = A(nodeListj, j);
              const Vector& gradAj = gradA(nodeListj, j);
              if (correctionOrder != CRKOrder::ZerothOrder) {
                Bj = B(nodeListj, j);
                gradBj = gradB(nodeListj, j);
              }
              if (correctionOrder == CRKOrder::QuadraticOrder) {
                Cj = C(nodeListj, j);
                gradCj = gradC(nodeListj, j);
              }
	      const DataType& Fj = fieldList(nodeListj, j);
	      GradientType& gradFj = result(nodeListj, j);

              // Node displacement.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;

              // Kernel weight and gradient.
              Scalar Wi, gWi, Wj, gWj;
              Vector gradWi, gradWj;
              CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, correctionOrder,  rij,  etai, Hi, Hdeti,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, -1e100, 1e100);
              CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, correctionOrder, -rij, -etaj, Hj, Hdetj, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj, -1e100, 1e100);

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
}

