set(CRKSPH_sources2d
    CRKSPHHydroBaseInst.cc.2.cc
    CRKSPHHydroBaseRZ.cc
    CRKSPHVariantInst.cc.2.cc
    ContinuityVolumePolicyInst.cc.2.cc
    HVolumePolicyInst.cc.2.cc
    SolidCRKSPHHydroBaseInst.cc.2.cc
    SolidCRKSPHHydroBaseRZ.cc
    computeCRKSPHCorrectionsInst.cc.2.cc
    computeCRKSPHEvaluationInst.cc.2.cc
    computeCRKSPHIntegralInst.cc.2.cc
    computeCRKSPHMomentsInst.cc.2.cc
    computeCRKSPHSumMassDensityInst.cc.2.cc
    computeCRKSPHSumVolumeInst.cc.2.cc
    computeHVolumesInst.cc.2.cc
    computeHullVolumesInst.cc.2.cc
    computeNeighborHullInst.cc.2.cc
    computeOccupancyVolumeInst.cc.2.cc
    computeSolidCRKSPHSumMassDensityInst.cc.2.cc
    detectSurfaceInst.cc.2.cc
    flagSurfaceNeighborsInst.cc.2.cc
    gradientCRKSPHInst.cc.2.cc
    interpolateCRKSPHInst.cc.2.cc
    )

set(CRKSPH_sources3d
    CRKSPHHydroBaseInst.cc.3.cc
    CRKSPHVariantInst.cc.3.cc
    ContinuityVolumePolicyInst.cc.3.cc
    HVolumePolicyInst.cc.3.cc
    SolidCRKSPHHydroBaseInst.cc.3.cc
    computeCRKSPHCorrectionsInst.cc.3.cc
    computeCRKSPHEvaluationInst.cc.3.cc
    computeCRKSPHIntegralInst.cc.3.cc
    computeCRKSPHMomentsInst.cc.3.cc
    computeCRKSPHSumMassDensityInst.cc.3.cc
    computeCRKSPHSumVolumeInst.cc.3.cc
    computeHVolumesInst.cc.3.cc
    computeHullVolumesInst.cc.3.cc
    computeNeighborHullInst.cc.3.cc
    computeOccupancyVolumeInst.cc.3.cc
    computeSolidCRKSPHSumMassDensityInst.cc.3.cc
    detectSurfaceInst.cc.3.cc
    flagSurfaceNeighborsInst.cc.3.cc
    gradientCRKSPHInst.cc.3.cc
    interpolateCRKSPHInst.cc.3.cc
    )

set(CRKSPH_sources
    centerOfMass.cc
    computeVoronoiVolume1d.cc
    computeVoronoiVolume2d.cc
    computeVoronoiVolume3d.cc
    )

set(CRKSPH_headers
    CRKSPHCorrectionParams.hh
    CRKSPHHydroBase.hh
    CRKSPHHydroBaseInline.hh
    CRKSPHHydroBaseRZ.hh
    CRKSPHSpecificThermalEnergyPolicy.hh
    CRKSPHUtilities.hh
    CRKSPHUtilitiesInline.hh
    CRKSPHVariant.hh
    ContinuityVolumePolicy.hh
    HVolumePolicy.hh
    HullVolumePolicy.hh
    SolidCRKSPHHydroBase.hh
    SolidCRKSPHHydroBaseInline.hh
    SolidCRKSPHHydroBaseRZ.hh
    SolidCRKSPHHydroBaseRZInline.hh
    SurfaceNodeCoupling.hh
    centerOfMass.hh
    computeCRKSPHCorrections.hh
    computeCRKSPHEvaluation.hh
    computeCRKSPHIntegral.hh
    computeCRKSPHMoments.hh
    computeCRKSPHSumMassDensity.hh
    computeCRKSPHSumVolume.hh
    computeHVolumes.hh
    computeHullSumMassDensity.hh
    computeHullVolumes.hh
    computeNeighborHull.hh
    computeOccupancyVolume.hh
    computeSolidCRKSPHSumMassDensity.hh
    computeVoronoiCentroids.hh
    computeVoronoiVolume.hh
    detectSurface.hh
    flagSurfaceNeighbors.hh
    gradientCRKSPH.hh
    interpolateCRKSPH.hh
    volumeSpacing.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_CRKSPH
                SOURCES     ${CRKSPH_sources} ${CRKSPH_sources2d} ${CRKSPH_sources3d}
                HEADERS     ${CRKSPH_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_CRKSPH
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CRKSPH_headers}
        DESTINATION include/CRKSPH
        )
