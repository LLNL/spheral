set(Hydro_inst
    EntropyPolicy
    GammaPolicy
    NonSymmetricSpecificThermalEnergyPolicy
    PositionPolicy
    PressurePolicy
    SecondMomentHourglassControl
    SoundSpeedPolicy
    SpecificFromTotalThermalEnergyPolicy
    SpecificThermalEnergyPolicy
    SumVoronoiMassDensityPolicy
    ThirdMomentHourglassControl
    VolumePolicy
    VoronoiMassDensityPolicy
   )


set(Hydro_sources
    HydroFieldNames.cc
    )

instantiate(Hydro_inst Hydro_sources)

set(Hydro_headers
    EntropyPolicy.hh
    GammaPolicy.hh
    HybridMassDensityPolicyInline.hh
    HydroFieldNames.hh
    HydroInline.hh
    MashCorrectionPolicy.hh
    MeshScaledMassDensityPolicy.hh
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
    TotalHydro.hh
    TotalHydroInline.hh
    VolumePolicy.hh
    VoronoiHourglassControl.hh
    VoronoiHourglassControlInline.hh
    VoronoiMassDensityPolicy.hh
    entropyWeightingFunction.hh
    )

INCLUDE_DIRECTORIES(. ../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Hydro
                SOURCES     ${Hydro_sources}
                HEADERS     ${Hydro_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Hydro
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Hydro_headers}
        DESTINATION include/Hydro
        )
