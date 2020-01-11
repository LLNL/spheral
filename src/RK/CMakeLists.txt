set(RK_inst
    ContinuityVolumePolicy
    HVolumePolicy
    RKCorrections
    RKUtilities
    ReproducingKernel
    ReproducingKernelMethods
    computeHVolumes
    computeHullVolumes
    computeOccupancyVolume
    computeRKSumVolume
    computeRKVolumes
    gradientRK
    hessianRK
    interpolateRK
   )


set(RK_sources
    HullVolumePolicy.cc
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
    ContinuityVolumePolicy.hh
    ContinuityVolumePolicyRZ.hh
    HVolumePolicy.hh
    HullVolumePolicy.hh
    RKCoefficients.hh
    RKCorrectionParams.hh
    RKCorrections.hh
    RKFieldNames.hh
    RKUtilities.hh
    RKUtilitiesInline.hh
    ReproducingKernel.hh
    ReproducingKernelInline.hh
    ReproducingKernelMethods.hh
    ReproducingKernelMethodsInline.hh
    computeHVolumes.hh
    computeHullVolumes.hh
    computeOccupancyVolume.hh
    computeRKSumVolume.hh
    computeRKVolumes.hh
    computeVoronoiVolume.hh
    gradientRK.hh
    hessianRK.hh
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
