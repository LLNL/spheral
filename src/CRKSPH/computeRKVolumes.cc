//------------------------------------------------------------------------------
// Compute the RK volumes
//------------------------------------------------------------------------------
#include "computeRKVolumes.hh"
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
                 const CRKVolumeType volumeType,
                 FieldList<Dimension, int>& surfacePoint,
                 FieldList<Dimension, typename Dimension::Vector>& deltaCentroid,
                 FieldList<Dimension, std::vector<typename Dimension::Vector>>& etaVoidPoints,
                 FieldList<Dimension, typename Dimension::FacetedVolume>& cells,
                 FieldList<Dimension, std::vector<CellFaceFlag>>& cellFaceFlags,
                 FieldList<Dimension, typename Dimension::Scalar>& volume) {
  if (volumeType == CRKVolumeType::CRKMassOverDensity) {
    volume.assignFields(mass/massDensity);
  }
  else if (volumeType == CRKVolumeType::CRKSumVolume) {
    computeCRKSPHSumVolume(connectivityMap, mW, position, mass, H, volume);
  }
  else if (volumeType == CRKVolumeType::CRKVoronoiVolume) {
    volume.assignFields(mass/massDensity);
    computeVoronoiVolume(position, H, connectivityMap, damage,
                         vector<typename Dimension::FacetedVolume>(),               // no boundaries
                         vector<vector<typename Dimension::FacetedVolume> >(),      // no holes
                         vector<Boundary<Dimension>*>(this->boundaryBegin(),        // boundaries
                                                      this->boundaryEnd()),
                         FieldList<Dimension, typename Dimension::Scalar>(),        // no weights
                         mSurfacePoint, volume, mDeltaCentroid, mEtaVoidPoints,    // return values
                         mCells,                                                    // return cells
                         mCellFaceFlags);                                           // node cell multimaterial faces
  }
  else if (volumeType == CRKVolumeType::CRKHullVolume) {
    computeHullVolumes(connectivityMap, mW.kernelExtent(), position, H, volume);
  }
  else if (volumeType == CRKVolumeType::HVolume) {
    const Scalar nPerh = volume.nodeListPtrs()[0]->nodesPerSmoothingScale();
    computeHVolumes(nPerh, H, volume);
  }
  else {
    VERIFY2(false, "Unknown CRK volume weighting.");
  }
}

} // end namespace Spheral
