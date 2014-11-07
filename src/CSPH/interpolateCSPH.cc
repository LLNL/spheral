//---------------------------------Spheral++------------------------------------
// Compute the CSPH interpolate.
//------------------------------------------------------------------------------
#include "interpolateCSPH.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "CSPHUtilities.hh"

namespace Spheral {
namespace CSPHSpace {

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
FieldSpace::FieldList<Dimension, DataType>
interpolateCSPH(const FieldSpace::FieldList<Dimension, DataType>& fieldList,
                const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                const bool coupleNodeLists,
                const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                const KernelSpace::TableKernel<Dimension>& W) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(A.size() == numNodeLists);
  REQUIRE(B.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Prepare the result.
  FieldList<Dimension, DataType> result;
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr != fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, DataType>("interpolate" + (*fieldItr)->name(), (*fieldItr)->nodeList()));
  }

  // Walk the FluidNodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const int firstGhostNodei = A[nodeListi]->nodeList().firstGhostNode();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Scalar wi = weight(nodeListi, i);
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar& Ai = A(nodeListi, i);
      const Vector& Bi = B(nodeListi, i);
      const DataType& Fi = fieldList(nodeListi, i);
      DataType& resulti = result(nodeListi, i);

      // Add our self-contribution.
      const Scalar W0 = W.kernelValue(0.0, Hdeti);
      resulti += wi*Fi*W0*Ai;

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Walk the neighbor nodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        if (coupleNodeLists or (nodeListi == nodeListj)) {
      
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

              // Only proceed if this node pair has not been calculated yet.
              if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                           nodeListj, j,
                                                           firstGhostNodej)) {

                // Get the state for node j.
                const Scalar wj = weight(nodeListj, j);
                const Vector& rj = position(nodeListj, j);
                const SymTensor& Hj = H(nodeListj, j);
                const Scalar Hdetj = Hj.Determinant();
                const Scalar& Aj = A(nodeListj, j);
                const Vector& Bj = B(nodeListj, j);
                const DataType& Fj = fieldList(nodeListj, j);
                DataType& resultj = result(nodeListj, j);

                // Node displacement.
                const Vector rij = ri - rj;
                const Vector etai = Hi*rij;
                const Vector etaj = Hj*rij;

                // Kernel weight.
                const Scalar Wj = CSPHKernel(W,  rij, etai, Hdeti, etaj, Hdetj, Ai, Bi);
                const Scalar Wi = CSPHKernel(W, -rij, etaj, Hdetj, etai, Hdeti, Aj, Bj);

                // Increment the pair-wise values.
                resulti += wj*Fj*Wj;
                resultj += wi*Fi*Wi;

              }
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

