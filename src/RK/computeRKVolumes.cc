//------------------------------------------------------------------------------
// Compute the RK volumes
//------------------------------------------------------------------------------
#include "computeRKVolumes.hh"

#include <vector>
#include "computeVoronoiVolume.hh"
#include "computeHullVolumes.hh"
#include "computeRKSumVolume.hh"
#include "computeHVolumes.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the volumes
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeRKVolumes(const ConnectivityMap<Dimension>& connectivityMap,
                 const TableKernel<Dimension>& W,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, typename Dimension::Scalar>& mass,
                 const FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                 const FieldList<Dimension, typename Dimension::SymTensor>& H,
                 const FieldList<Dimension, typename Dimension::SymTensor>& damage,
                 const std::vector<Boundary<Dimension>*>& boundaryConditions,
                 const RKVolumeType volumeType,
                 FieldList<Dimension, int>& surfacePoint,
                 FieldList<Dimension, typename Dimension::Vector>& deltaCentroid,
                 FieldList<Dimension, std::vector<typename Dimension::Vector>>& etaVoidPoints,
                 FieldList<Dimension, typename Dimension::FacetedVolume>& cells,
                 FieldList<Dimension, std::vector<CellFaceFlag>>& cellFaceFlags,
                 FieldList<Dimension, typename Dimension::Scalar>& volume) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  
  if (volumeType == RKVolumeType::RKMassOverDensity) {
    volume.assignFields(mass/massDensity);
  }
  else if (volumeType == RKVolumeType::RKSumVolume) {
    computeRKSumVolume(connectivityMap, W, position, mass, H, volume);
  }
  else if (volumeType == RKVolumeType::RKVoronoiVolume) {
    std::vector<std::vector<FacetedVolume>> holes;
    std::vector<FacetedVolume> facetedBoundaries;
    FieldList<Dimension, typename Dimension::Scalar> weights;
    volume.assignFields(mass/massDensity);
    computeVoronoiVolume(position, H, connectivityMap, damage,
                         facetedBoundaries, holes, boundaryConditions, weights,
                         surfacePoint, volume, deltaCentroid, etaVoidPoints,
                         cells, cellFaceFlags); 
  }
  else if (volumeType == RKVolumeType::RKHullVolume) {
    computeHullVolumes(connectivityMap, W.kernelExtent(), position, H, volume);
  }
  else if (volumeType == RKVolumeType::HVolume) {
    const Scalar nPerh = volume.nodeListPtrs()[0]->nodesPerSmoothingScale();
    computeHVolumes(nPerh, H, volume);
  }
  else {
    VERIFY2(false, "Unknown RK volume weighting.");
  }
}

} // end namespace Spheral
