include_directories(.)
set(Damage_inst
    DamageModel
    TensorDamageModel
    StrainPolicy
    TensorStrainPolicy
    TensorDamagePolicy
    YoungsModulusPolicy
    LongitudinalSoundSpeedPolicy
    weibullFlawDistribution
    computeFragmentField
    JohnsonCookDamage
    JohnsonCookFailureStrainPolicy
    JohnsonCookDamagePolicy
   )


set(Damage_sources )

instantiate(Damage_inst Damage_sources)

set(Damage_headers
    computeFragmentField.hh
    DamagedSoundSpeedPolicy.hh
    DamageGradientPolicy.hh
    DamageModel.hh
    DamageModelInline.hh
    EffectiveTensorDamagePolicy.hh
    GradyKippScalarDamage.hh
    JohnsonCookDamage.hh
    JohnsonCookDamageInline.hh
    JohnsonCookDamagePolicy.hh
    JohnsonCookFailureStrainPolicy.hh
    LongitudinalSoundSpeedPolicy.hh
    oneMinusDamage.hh
    ScalarDamageModel.hh
    ScalarDamagePolicy.hh
    StrainPolicy.hh
    TensorDamageModel.hh
    TensorDamageModelInline.hh
    TensorDamagePolicy.hh
    TensorStrainPolicy.hh
    weibullFlawDistribution.hh
    YoungsModulusPolicy.hh
    )

set(Damage_python_targets
    identifyFragments.py
    WeibullDamage.py
    GradyKippTensorDamage.py
    JohnsonCookDamageFactories.py
    )

blt_add_library(NAME        Spheral_Damage
                SOURCES     ${Damage_sources}
                HEADERS     ${Damage_headers}
                DEPENDS_ON  ${spheral_depends}
                SHARED      TRUE
                )

install(TARGETS              Spheral_Damage
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Damage_headers}
        DESTINATION include/Damage
        )

if(NOT ENABLE_CXXONLY)
  install(FILES       ${Damage_python_targets}
          DESTINATION .
          )
endif()

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Damage)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
