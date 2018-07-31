set(Hydro_sources2d
    EntropyPolicyInst.cc.2.cc
    GammaPolicyInst.cc.2.cc
    NonSymmetricSpecificThermalEnergyPolicyInst.cc.2.cc
    PositionPolicyInst.cc.2.cc
    PressurePolicyInst.cc.2.cc
    SecondMomentHourglassControlInst.cc.2.cc
    SoundSpeedPolicyInst.cc.2.cc
    SpecificFromTotalThermalEnergyPolicyInst.cc.2.cc
    SpecificThermalEnergyPolicyInst.cc.2.cc
    SumVoronoiMassDensityPolicyInst.cc.2.cc
    ThirdMomentHourglassControlInst.cc.2.cc
    VolumePolicyInst.cc.2.cc
    VoronoiMassDensityPolicyInst.cc.2.cc
    )

set(Hydro_sources3d
    EntropyPolicyInst.cc.3.cc
    GammaPolicyInst.cc.3.cc
    NonSymmetricSpecificThermalEnergyPolicyInst.cc.3.cc
    PositionPolicyInst.cc.3.cc
    PressurePolicyInst.cc.3.cc
    SecondMomentHourglassControlInst.cc.3.cc
    SoundSpeedPolicyInst.cc.3.cc
    SpecificFromTotalThermalEnergyPolicyInst.cc.3.cc
    SpecificThermalEnergyPolicyInst.cc.3.cc
    SumVoronoiMassDensityPolicyInst.cc.3.cc
    ThirdMomentHourglassControlInst.cc.3.cc
    VolumePolicyInst.cc.3.cc
    VoronoiMassDensityPolicyInst.cc.3.cc
    )

set(Hydro_sources
    HydroFieldNames.cc
    )

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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Hydro
                SOURCES     ${Hydro_sources2d} ${Hydro_sources3d} ${Hydro_sources}
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
