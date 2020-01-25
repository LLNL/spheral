set(RK_inst
    ContinuityVolumePolicy
    HVolumePolicy
    RKCorrections
    RKUtilities
    computeHVolumes
    computeHullVolumes
    computeOccupancyVolume
    computeRKSumVolume
    computeRKVolumes
    interpolateRK
   )


set(RK_sources
    HullVolumePolicy.cc
    computeVoronoiVolume.cc
    computeVoronoiVolume1d.cc
   )

if(ENABLE_2D)
  list(APPEND RK_sources
       ContinuityVolumePolicyRZ.cc
      )
endif()

instantiate(RK_inst RK_sources)

set(RK_headers
    ContinuityVolumePolicy.hh
    ContinuityVolumePolicyRZ.hh
    HVolumePolicy.hh
    HullVolumePolicy.hh
    RKCorrectionParams.hh
    RKCorrections.hh
    RKUtilities.hh
    RKUtilitiesInline.hh
    computeHVolumes.hh
    computeHullVolumes.hh
    computeOccupancyVolume.hh
    computeRKSumVolume.hh
    computeRKVolumes.hh
    computeVoronoiVolume.hh
    interpolateRK.hh
    )

blt_add_library(NAME        Spheral_RK
                SOURCES     ${RK_sources}
                HEADERS     ${RK_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_RK
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${RK_headers}
        DESTINATION include/RK
        )
