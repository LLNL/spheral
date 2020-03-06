include_directories(.)
set(RK_inst
  computeRKSumVolume
  computeOccupancyVolume
  computeHullVolumes
  computeHVolumes
  HVolumePolicy
  ContinuityVolumePolicy
  computeRKVolumes
  RKCorrections
  RKUtilities
  ReproducingKernelMethods
  ReproducingKernel
  interpolateRK
  gradientRK
  hessianRK
  )

set(RK_sources
  computeVoronoiVolume.cc
  computeVoronoiVolume1d.cc
  RKFieldNames.cc
  )

if(ENABLE_2D)
  list(APPEND RK_sources
    ContinuityVolumePolicyRZ.cc
    )
endif()

instantiate(RK_inst RK_sources)

set(RK_headers
  computeHullVolumes.hh
  computeHVolumes.hh
  computeOccupancyVolume.hh
  computeRKSumVolume.hh
  computeRKVolumes.hh
  computeVoronoiVolume.hh
  ContinuityVolumePolicy.hh
  ContinuityVolumePolicyRZ.hh
  HullVolumePolicy.hh
  HVolumePolicy.hh
  interpolateRK.hh
  RKCorrectionParams.hh
  RKCorrections.hh
  RKUtilities.hh
  RKUtilitiesInline.hh
  )

spheral_add_cxx_library(RK)

install(TARGETS              Spheral_RK
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${RK_headers}
        DESTINATION include/RK
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_RK)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
