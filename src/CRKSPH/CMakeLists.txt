set(CRKSPH_inst
    CRKSPHHydroBase
    ContinuityVolumePolicy
    SolidCRKSPHHydroBase
    computeCRKSPHCorrections
    computeCRKSPHEvaluation
    computeCRKSPHIntegral
    computeCRKSPHMoments
    computeCRKSPHSumMassDensity
    computeCRKSPHSumVolume
    computeHVolumes
    computeHullVolumes
    computeRKVolumes
    computeSolidCRKSPHSumMassDensity
    detectSurface
    editMultimaterialSurfaceTopology
    gradientCRKSPH
    interpolateCRKSPH
    zerothOrderSurfaceCorrections
    RKCorrections
    RKUtilities
   )


set(CRKSPH_sources
    centerOfMass.cc
    ContinuityVolumePolicyRZ.cc
    computeVoronoiVolume.cc
    computeVoronoiVolume1d.cc
   )

if(ENABLE_2D)
  list(APPEND CRKSPH_sources
       CRKSPHHydroBaseRZ.cc
       SolidCRKSPHHydroBaseRZ.cc
      )
endif()

instantiate(CRKSPH_inst CRKSPH_sources)

set(CRKSPH_headers
    CRKSPHCorrectionParams.hh
    CRKSPHHydroBase.hh
    CRKSPHHydroBaseInline.hh
    CRKSPHHydroBaseRZ.hh
    CRKSPHUtilities.hh
    CRKSPHUtilitiesInline.hh
    ContinuityVolumePolicy.hh
    ContinuityVolumePolicyRZ.hh
    HullVolumePolicy.hh
    SolidCRKSPHHydroBase.hh
    SolidCRKSPHHydroBaseInline.hh
    SolidCRKSPHHydroBaseRZ.hh
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
    computeRKVolumes.hh
    computeSolidCRKSPHSumMassDensity.hh
    computeVoronoiCentroids.hh
    computeVoronoiVolume.hh
    detectSurface.hh
    editMultimaterialSurfaceTopology.hh
    gradientCRKSPH.hh
    interpolateCRKSPH.hh
    volumeSpacing.hh
    zerothOrderSurfaceCorrections.hh
    RKCorrections.hh
    RKUtilities.hh
    RKUtilitiesInline.hh
    )

blt_add_library(NAME        Spheral_CRKSPH
                SOURCES     ${CRKSPH_sources}
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
