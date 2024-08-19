include_directories(.)

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

# Add explicit instantiations
list(APPEND RK_sources
     computeRKSumVolume
     computeOccupancyVolume
     computeHullVolumes
     computeHVolumes
     HVolumePolicy
     ContinuityVolumePolicy
     computeRKVolumes
     RKCorrections
     RKUtilities012
     RKUtilities3
     RKUtilities4
     RKUtilities5
     RKUtilities6
     RKUtilities7
     ReproducingKernelMethods
     ReproducingKernel
     interpolateRK
     gradientRK
     hessianRK
     )

set(RK_headers
    ContinuityVolumePolicy.hh
    ContinuityVolumePolicyRZ.hh
    HVolumePolicy.hh
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

spheral_add_obj_library(RK SPHERAL_OBJ_LIBS)

