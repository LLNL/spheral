#include "computeCRKSPHSumVolume.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Function to compute the per dimension volume multiplier.
//------------------------------------------------------------------------------
namespace {
template<typename Dimension> double volumeElement();

template<> double volumeElement<Dim<1> >() {
  return 2.0;
}
  
template<> double volumeElement<Dim<2> >() {
  return M_PI;
}
  
template<> double volumeElement<Dim<3> >() {
  return 4.0/3.0*M_PI;
}
}

//------------------------------------------------------------------------------
// Compute the CRKSPH volume summation.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeCRKSPHSumVolume(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::Scalar>& mass,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Scalar>& vol) {

  // Pre-conditions.
  const size_t numNodeLists = vol.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // Get the maximum allowed volume in eta space.
  const auto etaVolMax = Dimension::pownu(0.5*W.kernelExtent()) * volumeElement<Dimension>();

  // Zero it out.
  vol = 0.0;

  // For our first pass compute the effective volume per point.
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = connectivityMap.numNodes(nodeListi);

    // Iterate over the nodes in this node list.
    for (auto iItr = 0; iItr < ni; ++iItr) {
      const auto i = connectivityMap.ithNode(nodeListi, iItr);

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const auto& connectivity = fullConnectivity[nodeListj];
        for (auto jItr = connectivity.begin(); jItr != connectivity.end(); ++jItr) {
          const auto  j = *jItr;
          const auto& rj = position(nodeListj, j);
          const auto& Hj = H(nodeListj, j);
          const auto  Hdetj = Hj.Determinant();
          vol(nodeListi, i) += W.kernelValue((Hi*(ri - rj)).magnitude(), Hdeti);
        }
      }
  
      // Add the self-contribution.
      vol(nodeListi, i) += W.kernelValue(0.0, Hdeti);
      CHECK(vol(nodeListi, i) > 0.0);
      vol(nodeListi, i) = min(etaVolMax/Hdeti, 1.0/vol(nodeListi, i));
    }
  }
}

}

