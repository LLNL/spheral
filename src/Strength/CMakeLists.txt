set(Strength_inst
    BulkModulusPolicy
    DeviatoricStressPolicy
    MeltEnergyPolicy
    PlasticStrainPolicy
    ShearModulusPolicy
    StrengthSoundSpeedPolicy
    YieldStrengthPolicy
   )

set(Strength_sources
    SolidFieldNames.cc
    )

if(ENABLE_2D)
  list(APPEND Strength_sources
       RZPlasticStrainPolicy.cc
      )
endif()

instantiate(Strength_inst Strength_sources)

set(Strength_headers
    BulkModulusPolicy.hh
    DeviatoricStressPolicy.hh
    MeltEnergyPolicy.hh
    PlasticStrainPolicy.hh
    RZPlasticStrainPolicy.hh
    ShearModulusPolicy.hh
    SolidFieldNames.hh
    StrengthSoundSpeedPolicy.hh
    YieldStrengthPolicy.hh
    effectiveKernelVolume.hh
    )

blt_add_library(NAME        Spheral_Strength
                SOURCES     ${Strength_sources}
                HEADERS     ${Strength_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Strength
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Strength_headers}
        DESTINATION include/Strength
        )
