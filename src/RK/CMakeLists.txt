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
     computeRKSumVolumeInst.cc
     computeOccupancyVolumeInst.cc
     computeHullVolumesInst.cc
     computeHVolumesInst.cc
     HVolumePolicyInst.cc
     ContinuityVolumePolicyInst.cc
     computeRKVolumesInst.cc
     RKCorrectionsInst.cc
     RKUtilities012Inst.cc
     RKUtilities3Inst.cc
     RKUtilities4Inst.cc
     RKUtilities5Inst.cc
     RKUtilities6Inst.cc
     RKUtilities7Inst.cc
     ReproducingKernelMethodsInst.cc
     ReproducingKernelInst.cc
     interpolateRKInst.cc
     gradientRKInst.cc
     hessianRKInst.cc
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

