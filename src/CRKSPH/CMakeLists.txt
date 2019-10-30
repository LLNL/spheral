set(CRKSPH_inst
  computeCRKSPHCorrections
  computeCRKSPHEvaluation
  computeCRKSPHIntegral
  computeCRKSPHMoments
  computeCRKSPHSumMassDensity
  computeCRKSPHSumVolume
  #computeHullSumMassDensity
  computeHullVolumes
  computeHVolumes
  computeNeighborHull
  computeOccupancyVolume
  computeSolidCRKSPHSumMassDensity
  ContinuityVolumePolicy
  #CRKSPHHydroBase
  #CRKSPHVariant
  detectSurface
  editMultimaterialSurfaceTopology
  gradientCRKSPH
  HVolumePolicy
  interpolateCRKSPH
  #SolidCRKSPHHydroBase
  zerothOrderSurfaceCorrections
  )

set(CRKSPH_sources 
  centerOfMass.cc
  computeVoronoiVolume1d.cc
  computeVoronoiVolume2d.cc
  computeVoronoiVolume3d.cc
  )
instantiate(CRKSPH_inst CRKSPH_sources)

set(CRKSPH_headers
  centerOfMass.hh
  computeCRKSPHCorrections.hh
  computeCRKSPHEvaluation.hh
  computeCRKSPHIntegral.hh
  computeCRKSPHMoments.hh
  computeCRKSPHSumMassDensity.hh
  computeCRKSPHSumVolume.hh
  #computeHullSumMassDensity.hh
  computeHullVolumes.hh
  computeHVolumes.hh
  computeNeighborHull.hh
  computeOccupancyVolume.hh
  computeSolidCRKSPHSumMassDensity.hh
  computeVoronoiCentroids.hh
  computeVoronoiVolume.hh
  ContinuityVolumePolicy.hh
  CRKSPHCorrectionParams.hh
  CRKSPHHydroBase.hh
  CRKSPHHydroBaseInline.hh
  CRKSPHHydroBaseRZ.hh
  CRKSPHSpecificThermalEnergyPolicy.hh
  CRKSPHUtilities.hh
  CRKSPHUtilitiesInline.hh
  CRKSPHVariant.hh
  detectSurface.hh
  editMultimaterialSurfaceTopology.hh
  gradientCRKSPH.hh
  HullVolumePolicy.hh
  HVolumePolicy.hh
  interpolateCRKSPH.hh
  SolidCRKSPHHydroBase.hh
  SolidCRKSPHHydroBaseInline.hh
  SolidCRKSPHHydroBaseRZ.hh
  SolidCRKSPHHydroBaseRZInline.hh
  SurfaceNodeCoupling.hh
  volumeSpacing.hh
  zerothOrderSurfaceCorrections.hh
  )


set(CRKSPH_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/CRKSPH
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_CRKSPH
                HEADERS     ${CRKSPH_headers}
                SOURCES     ${CRKSPH_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${CRKSPH_includes}
                )

install(TARGETS             Spheral_CRKSPH
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CRKSPH_headers}
        DESTINATION include/CRKSPH
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_CRKSPH)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
