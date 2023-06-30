include_directories(.)
set(Hydro_inst
    EntropyPolicy
    CompatibleDifferenceSpecificThermalEnergyPolicy
    PositionPolicy
    PressurePolicy
    SoundSpeedPolicy
    SpecificFromTotalThermalEnergyPolicy
    SpecificThermalEnergyPolicy
    SumVoronoiMassDensityPolicy
    VolumePolicy
    VoronoiMassDensityPolicy
    GammaPolicy
    SecondMomentHourglassControl
    ThirdMomentHourglassControl
   )


set(Hydro_sources
    HydroFieldNames.cc
    )

if (ENABLE_1D AND ENABLE_2D)
  list(APPEND Hydro_sources SphericalPositionPolicy.cc)
endif()

instantiate(Hydro_inst Hydro_sources)

set(Hydro_headers
    EntropyPolicy.hh
    GammaPolicy.hh
    HybridMassDensityPolicyInline.hh
    HydroFieldNames.hh
    HydroInline.hh
    MashCorrectionPolicy.hh
    MeshScaledMassDensityPolicy.hh
    CompatibleDifferenceSpecificThermalEnergyPolicy.hh
    NonSymmetricSpecificThermalEnergyPolicy.hh
    PositionPolicy.hh
    PressurePolicy.hh
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

spheral_add_obj_library(Hydro)

