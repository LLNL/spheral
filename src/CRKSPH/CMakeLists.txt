include_directories(.)
set(CRKSPH_inst
    computeCRKSPHSumMassDensity
    #computeCRKSPHSumVolume
    #computeOccupancyVolume
    computeSolidCRKSPHSumMassDensity
    computeCRKSPHMoments
    computeCRKSPHCorrections
    #computeCRKSPHEvaluation
    #computeCRKSPHIntegral
    #computeHullVolumes
    #computeHVolumes
    #computeNeighborHull
    CRKSPHHydroBase
    #CRKSPHVariant
    interpolateCRKSPH
    gradientCRKSPH
    #HVolumePolicy
    detectSurface
    #ContinuityVolumePolicy
    SolidCRKSPHHydroBase
    editMultimaterialSurfaceTopology
    zerothOrderSurfaceCorrections
   )


set(CRKSPH_sources
    centerOfMass.cc
    #computeVoronoiVolume.cc
    #computeVoronoiVolume1d.cc
   )

if(ENABLE_2D)
  list(APPEND CRKSPH_sources
       CRKSPHHydroBaseRZ.cc
       #computeVoronoiVolume2d.cc
       SolidCRKSPHHydroBaseRZ.cc
       #ContinuityVolumePolicyRZ.cc
      )
endif()

if(ENABLE_3D)
  list(APPEND CRKSPH_sources
      #computeVoronoiVolume3d.cc
      )
endif()

instantiate(CRKSPH_inst CRKSPH_sources)

set(CRKSPH_headers
  #CRKSPHCorrectionParams.hh
    CRKSPHHydroBase.hh
    CRKSPHHydroBaseInline.hh
    CRKSPHHydroBaseRZ.hh
    #CRKSPHSpecificThermalEnergyPolicy.hh
    CRKSPHUtilities.hh
    CRKSPHUtilitiesInline.hh
    #ContinuityVolumePolicy.hh
    #HullVolumePolicy.hh
    SolidCRKSPHHydroBase.hh
    SolidCRKSPHHydroBaseInline.hh
    SolidCRKSPHHydroBaseRZ.hh
    #SolidCRKSPHHydroBaseRZInline.hh
    SurfaceNodeCoupling.hh
    centerOfMass.hh
    computeCRKSPHCorrections.hh
    computeCRKSPHEvaluation.hh
    computeCRKSPHIntegral.hh
    computeCRKSPHMoments.hh
    computeCRKSPHSumMassDensity.hh
    #computeCRKSPHSumVolume.hh
    #computeHVolumes.hh
    computeHullSumMassDensity.hh
    #computeHullVolumes.hh
    computeSolidCRKSPHSumMassDensity.hh
    computeVoronoiCentroids.hh
    #computeVoronoiVolume.hh
    detectSurface.hh
    #flagSurfaceNeighbors.hh
    gradientCRKSPH.hh
    interpolateCRKSPH.hh
    volumeSpacing.hh
    zerothOrderSurfaceCorrections.hh
    )

set(CRKSPH_python_targets
  CRKSPHHydros.py
  )

blt_add_library(NAME        Spheral_CRKSPH
                SOURCES     ${CRKSPH_sources}
                HEADERS     ${CRKSPH_headers}
                DEPENDS_ON  ${spheral_depends}
                SHARED      TRUE
                )

install(TARGETS              Spheral_CRKSPH
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CRKSPH_headers}
        DESTINATION include/CRKSPH
        )

install(FILES       ${CRKSPH_python_targets}
        DESTINATION .
        )
    
list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_CRKSPH)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
