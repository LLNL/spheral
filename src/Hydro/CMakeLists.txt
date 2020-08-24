include_directories(.)
set(Hydro_inst
    EntropyPolicy
    NonSymmetricSpecificThermalEnergyPolicy
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

if (ENABLE_2D)
  list(APPEND Hydro_sources RZNonSymmetricSpecificThermalEnergyPolicy.cc)
endif()

instantiate(Hydro_inst Hydro_sources)

set(Hydro_headers
    EntropyPolicy.hh
    HybridMassDensityPolicyInline.hh
    HydroFieldNames.hh
    HydroInline.hh
    MashCorrectionPolicy.hh
    MeshScaledMassDensityPolicy.hh
    NonSymmetricSpecificThermalEnergyPolicy.hh
    PositionPolicy.hh
    PressurePolicy.hh
    SoundSpeedPolicy.hh
    SpecificFromTotalThermalEnergyPolicy.hh
    SpecificThermalEnergyPolicy.hh
    SumVoronoiMassDensityPolicy.hh
    VolumePolicy.hh
    VoronoiMassDensityPolicy.hh
    entropyWeightingFunction.hh
    )

spheral_add_cxx_library(Hydro)

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Hydro)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
