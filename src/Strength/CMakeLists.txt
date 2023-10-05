include_directories(.)
set(Strength_inst
    BulkModulusPolicy
    DeviatoricStressPolicy
    MeltEnergyPolicy
    PlasticStrainPolicy
    ShearModulusPolicy
    StrengthSoundSpeedPolicy
    YieldStrengthPolicy
    PorousPressurePolicy
    PorousSoundSpeedPolicy
    PorousGammaPolicy
    PorousBulkModulusPolicy
    PorousEntropyPolicy
    PorousShearModulusPolicy
    PorousYieldStrengthPolicy
   )

set(Strength_sources
    SolidFieldNames.cc
    )

instantiate(Strength_inst Strength_sources)

set(Strength_headers
    BulkModulusPolicy.hh
    DeviatoricStressPolicy.hh
    MeltEnergyPolicy.hh
    PlasticStrainPolicy.hh
    ShearModulusPolicy.hh
    SolidFieldNames.hh
    StrengthSoundSpeedPolicy.hh
    YieldStrengthPolicy.hh
    effectiveKernelVolume.hh
    PorousPressurePolicy.hh
    PorousSoundSpeedPolicy.hh
    PorousGammaPolicy.hh
    PorousBulkModulusPolicy.hh
    PorousEntropyPolicy.hh
    PorousShearModulusPolicy.hh
    PorousYieldStrengthPolicy.hh
    )

spheral_add_obj_library(Strength)

