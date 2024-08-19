//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "RK/computeRKVolumes.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

template void computeRKVolumes<Dim<1>>(const ConnectivityMap<Dim<1>>& connectivityMap,
                                              const TableKernel<Dim<1>>& W,
                                              const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                              const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
                                              const FieldList<Dim<1>, Dim<1>::Scalar>& massDensity,
                                              const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                                              const FieldList<Dim<1>, Dim<1>::SymTensor>& damage,
                                              const std::vector<Dim<1>::FacetedVolume>& facetedBoundaries,
                                              const std::vector<std::vector<Dim<1>::FacetedVolume>>& facetedHoles,
                                              const std::vector<Boundary<Dim<1>>*>& boundaryConditions,
                                              const RKVolumeType volumeType,
                                              FieldList<Dim<1>, int>& surfacePoint,
                                              FieldList<Dim<1>, Dim<1>::Vector>& deltaCentroid,
                                              FieldList<Dim<1>, std::vector<Dim<1>::Vector>>& etaVoidPoints,
                                              FieldList<Dim<1>, Dim<1>::FacetedVolume>& cells,
                                              FieldList<Dim<1>, std::vector<CellFaceFlag>>& cellFaceFlags,
                                              FieldList<Dim<1>, Dim<1>::Scalar>& volume);

#endif

#if defined(SPHERAL_ENABLE_2D)

template void computeRKVolumes<Dim<2>>(const ConnectivityMap<Dim<2>>& connectivityMap,
                                              const TableKernel<Dim<2>>& W,
                                              const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                              const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
                                              const FieldList<Dim<2>, Dim<2>::Scalar>& massDensity,
                                              const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                                              const FieldList<Dim<2>, Dim<2>::SymTensor>& damage,
                                              const std::vector<Dim<2>::FacetedVolume>& facetedBoundaries,
                                              const std::vector<std::vector<Dim<2>::FacetedVolume>>& facetedHoles,
                                              const std::vector<Boundary<Dim<2>>*>& boundaryConditions,
                                              const RKVolumeType volumeType,
                                              FieldList<Dim<2>, int>& surfacePoint,
                                              FieldList<Dim<2>, Dim<2>::Vector>& deltaCentroid,
                                              FieldList<Dim<2>, std::vector<Dim<2>::Vector>>& etaVoidPoints,
                                              FieldList<Dim<2>, Dim<2>::FacetedVolume>& cells,
                                              FieldList<Dim<2>, std::vector<CellFaceFlag>>& cellFaceFlags,
                                              FieldList<Dim<2>, Dim<2>::Scalar>& volume);

#endif

#if defined(SPHERAL_ENABLE_3D)

template void computeRKVolumes<Dim<3>>(const ConnectivityMap<Dim<3>>& connectivityMap,
                                              const TableKernel<Dim<3>>& W,
                                              const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                              const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                              const FieldList<Dim<3>, Dim<3>::Scalar>& massDensity,
                                              const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                                              const FieldList<Dim<3>, Dim<3>::SymTensor>& damage,
                                              const std::vector<Dim<3>::FacetedVolume>& facetedBoundaries,
                                              const std::vector<std::vector<Dim<3>::FacetedVolume>>& facetedHoles,
                                              const std::vector<Boundary<Dim<3>>*>& boundaryConditions,
                                              const RKVolumeType volumeType,
                                              FieldList<Dim<3>, int>& surfacePoint,
                                              FieldList<Dim<3>, Dim<3>::Vector>& deltaCentroid,
                                              FieldList<Dim<3>, std::vector<Dim<3>::Vector>>& etaVoidPoints,
                                              FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells,
                                              FieldList<Dim<3>, std::vector<CellFaceFlag>>& cellFaceFlags,
                                              FieldList<Dim<3>, Dim<3>::Scalar>& volume);

#endif
}