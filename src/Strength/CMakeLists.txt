include_directories(.)
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
       #RZDeviatoricStressPolicy.cc
       RZPlasticStrainPolicy.cc
      )
endif()

instantiate(Strength_inst Strength_sources)

set(Strength_headers
    BulkModulusPolicy.hh
    DeviatoricStressPolicy.hh
    MeltEnergyPolicy.hh
    PlasticStrainPolicy.hh
    #RZDeviatoricStressPolicy.hh
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
                SHARED      TRUE
                )

install(TARGETS              Spheral_Strength
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Strength_headers}
        DESTINATION include/Strength
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Strength)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
