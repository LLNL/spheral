//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH interpolate.
//------------------------------------------------------------------------------
#include "interpolateCRKSPH.hh"
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
FieldList<Dimension, DataType>
interpolateCRKSPH(const FieldList<Dimension, DataType>& fieldList,
                  const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::Scalar>& weight,
                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                  const FieldList<Dimension, typename Dimension::Scalar>& A,
                  const FieldList<Dimension, typename Dimension::Vector>& B,
                  const FieldList<Dimension, typename Dimension::Tensor>& C,
                  const ConnectivityMap<Dimension>& connectivityMap,
                  const CRKOrder correctionOrder,
                  const TableKernel<Dimension>& W,
                  const NodeCoupling& nodeCoupling) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(A.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists or correctionOrder == CRKOrder::ZerothOrder);
  REQUIRE(C.size() == numNodeLists or correctionOrder != CRKOrder::QuadraticOrder);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Prepare the result.
  FieldList<Dimension, DataType> result;
  result.copyFields();
  for (auto fieldItr = fieldList.begin();
       fieldItr != fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DataType>("interpolate" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
  }

  // Walk the FluidNodeLists.
  auto Bi = Vector::zero, Bj = Vector::zero;
  auto Ci = Tensor::zero, Cj = Tensor::zero;
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
      if (correctionOrder != CRKOrder::ZerothOrder) Bi = B(nodeListi, i);
      if (correctionOrder == CRKOrder::QuadraticOrder) Ci = C(nodeListi, i);
      const auto& Fi = fieldList(nodeListi, i);
      auto& resulti = result(nodeListi, i);

      // Add our self-contribution.
      const auto W0 = W.kernelValue(0.0, Hdeti);
      resulti += weight(nodeListi, i)*Fi*W0*Ai;

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
              const auto wi = fij*weight(nodeListi, i);
              const auto wj = fij*weight(nodeListj, j);

              // Get the state for node j.
              const auto& rj = position(nodeListj, j);
              const auto& Hj = H(nodeListj, j);
              const auto  Hdetj = Hj.Determinant();
              const auto  Aj = A(nodeListj, j);
              if (correctionOrder != CRKOrder::ZerothOrder) Bj = B(nodeListj, j);
              if (correctionOrder == CRKOrder::QuadraticOrder) Cj = C(nodeListj, j);
              const auto& Fj = fieldList(nodeListj, j);
              auto& resultj = result(nodeListj, j);

              // Node displacement.
              const auto rij = ri - rj;
              const auto etai = Hi*rij;
              const auto etaj = Hj*rij;

              // Kernel weight.
              const auto Wj = CRKSPHKernel(W, correctionOrder,  rij,  etai, Hdeti,  etaj, Hdetj, Ai, Bi, Ci);
              const auto Wi = CRKSPHKernel(W, correctionOrder, -rij, -etaj, Hdetj, -etai, Hdeti, Aj, Bj, Cj);

              // Increment the pair-wise values.
              resulti += wj*Fj*Wj;
              resultj += wi*Fi*Wi;

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
