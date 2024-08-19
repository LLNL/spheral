include_directories(.)

set(Hydro_sources
    HydroFieldNames.cc
    )

if (ENABLE_2D)
  list(APPEND Hydro_sources RZNonSymmetricSpecificThermalEnergyPolicy.cc)
endif()

if (ENABLE_1D AND ENABLE_2D)
  list(APPEND Hydro_sources SphericalPositionPolicy.cc)
endif()

# Add explicit instantiations
list(APPEND Hydro_sources
     EntropyPolicyInst.cc
     CompatibleDifferenceSpecificThermalEnergyPolicyInst.cc
     NonSymmetricSpecificThermalEnergyPolicyInst.cc
     PressurePolicyInst.cc
     SoundSpeedPolicyInst.cc
     SpecificFromTotalThermalEnergyPolicyInst.cc
     SpecificThermalEnergyPolicyInst.cc
     SumVoronoiMassDensityPolicyInst.cc
     VolumePolicyInst.cc
     VoronoiMassDensityPolicyInst.cc
     GammaPolicyInst.cc
     SecondMomentHourglassControlInst.cc
     ThirdMomentHourglassControlInst.cc
     )

set(Hydro_headers
    EntropyPolicy.hh
    GammaPolicy.hh
    HydroFieldNames.hh
    HydroInline.hh
    CompatibleDifferenceSpecificThermalEnergyPolicy.hh
    NonSymmetricSpecificThermalEnergyPolicy.hh
    PressurePolicy.hh
    RZNonSymmetricSpecificThermalEnergyPolicy.hh
    SecondMomentHourglassControl.hh
    SecondMomentHourglassControlInline.hh
    SoundSpeedPolicy.hh
    SpecificFromTotalThermalEnergyPolicy.hh
    SpecificThermalEnergyPolicy.hh
    SumVoronoiMassDensityPolicy.hh
    ThirdMomentHourglassControl.hh
    ThirdMomentHourglassControlInline.hh
    VolumePolicy.hh
    VoronoiHourglassControl.hh
    VoronoiHourglassControlInline.hh
    VoronoiMassDensityPolicy.hh
    entropyWeightingFunction.hh
    )

spheral_add_obj_library(Hydro SPHERAL_OBJ_LIBS)

