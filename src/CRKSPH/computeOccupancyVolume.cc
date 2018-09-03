//------------------------------------------------------------------------------
// Compute the CRKSPH volume summation.
//------------------------------------------------------------------------------
#include "computeOccupancyVolume.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

template<typename Dimension>
void
computeOccupancyVolume(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Scalar>& vol) {

  // Pre-conditions.
  const size_t numNodeLists = vol.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  // Zero it out.
  vol = 0.0;

  // Extent of the kernel.
  const Scalar kernelExtent = W.kernelExtent();
  const Scalar volFactor = 2.0*kernelExtent;

  // Walk the NodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

    // Iterate over the nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const unsigned Ni = connectivityMap.numNeighborsForNode(nodeListi, i) + 1;
      vol(nodeListi, i) = volFactor/(Ni*Hdeti);
    }
  }
}

}

