//------------------------------------------------------------------------------
// Compute the Shepards interpolation of a FieldList.
//------------------------------------------------------------------------------
#include "computeShepardsInterpolation.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/DataTypeTraits.hh"

namespace Spheral {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;

template<typename Dimension, typename DataType>
FieldList<Dimension, DataType>
computeShepardsInterpolation(const FieldList<Dimension, DataType>& fieldList,
                             const ConnectivityMap<Dimension>& connectivityMap,
                             const TableKernel<Dimension>& W,
                             const FieldList<Dimension, typename Dimension::Vector>& position,
                             const FieldList<Dimension, typename Dimension::SymTensor>& H,
                             const FieldList<Dimension, typename Dimension::Scalar>& weight) {

  // Pre-conditions.
  const size_t numNodeLists = fieldList.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(weight.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  // Some scratch variables.
  size_t nodeListi, nodeListj;
  Scalar Wi, Wj;
  Vector rij, etai, etaj;
  Vector Bi = Vector::zero, Bj = Vector::zero;
  Tensor Ci = Tensor::zero, Cj = Tensor::zero;

  // Prepare the return value.
  FieldList<Dimension, DataType> result(FieldSpace::FieldStorageType::CopyFields);
  FieldList<Dimension, Scalar> wsum(FieldSpace::FieldStorageType::CopyFields);
  for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    result.appendNewField(fieldList[nodeListi]->name(), fieldList[nodeListi]->nodeList(), DataTypeTraits<DataType>::zero());
    wsum.appendNewField("weight sum", fieldList[nodeListi]->nodeList(), 0.0);
  }

  // Walk all the points.
  for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto  weighti = weight(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // Walk the neighbor NodeLists.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      for (nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const auto& connectivity = fullConnectivity[nodeListj];
        const int firstGhostNodej = fieldList[nodeListi]->nodeList().firstGhostNode();
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const int j = *jItr;

          // Check if this node pair has already been calculated.
          if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                       nodeListj, j,
                                                       firstGhostNodej)) {
            const auto& rj = position(nodeListi, j);
            const auto  weightj = weight(nodeListi, j);
            const auto& Hj = H(nodeListi, j);
            const auto  Hdetj = Hj.Determinant();

            // Kernel weighting and gradient.
            rij = ri - rj;
            etai = Hi*rij;
            etaj = Hj*rij;
            Wi = W.kernelValue(etai.magnitude(), Hdeti);
            Wj = W.kernelValue(etaj.magnitude(), Hdetj);

            // Sum the pair-wise contributions.
            wsum(nodeListi, i) += weightj*Wj;
            wsum(nodeListj, j) += weighti*Wi;
            result(nodeListi, i) += fieldList(nodeListi, i) * weightj*Wj;
            result(nodeListj, j) += fieldList(nodeListj, j) * weighti*Wi;
          }
        }
      }
  
      // Finalize the interpolation for node i.
      const auto wW0 = weighti*W.kernelValue(0.0, Hdeti);
      result(nodeListi, i) = ((result(nodeListi, i) + fieldList(nodeListi, i) * wW0)/
                              (wsum(nodeListi, i) + wW0));
    }
  }

  // That's it.
  return result;
}

}

