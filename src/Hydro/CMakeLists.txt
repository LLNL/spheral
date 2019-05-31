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
   )


set(Hydro_sources
    HydroFieldNames.cc
    )

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
    TotalHydro.hh
    TotalHydroInline.hh
    VolumePolicy.hh
    VoronoiMassDensityPolicy.hh
    entropyWeightingFunction.hh
    )

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
