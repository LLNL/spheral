set(Damage_inst
    DamageModel
    JohnsonCookDamage
    JohnsonCookDamagePolicy
    JohnsonCookFailureStrainPolicy
    LongitudinalSoundSpeedPolicy
    StrainPolicy
    TensorDamageModel
    TensorDamagePolicy
    TensorStrainPolicy
    YoungsModulusPolicy
    computeFragmentField
    weibullFlawDistribution
   )


set(Damage_sources )

instantiate(Damage_inst Damage_sources)

set(Damage_headers
    DamageGradientPolicy.hh
    DamageModel.hh
    DamageModelInline.hh
    DamagedSoundSpeedPolicy.hh
    EffectiveTensorDamagePolicy.hh
    GradyKippScalarDamage.hh
    JohnsonCookDamage.hh
    JohnsonCookDamageInline.hh
    JohnsonCookDamagePolicy.hh
    JohnsonCookFailureStrainPolicy.hh
    LongitudinalSoundSpeedPolicy.hh
    ScalarDamageModel.hh
    ScalarDamagePolicy.hh
    StrainPolicy.hh
    TensorDamageModel.hh
    TensorDamageModelInline.hh
    TensorDamagePolicy.hh
    TensorStrainPolicy.hh
    YoungsModulusPolicy.hh
    computeFragmentField.hh
    oneMinusDamage.hh
    weibullFlawDistribution.hh
    )

INCLUDE_DIRECTORIES(. ../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Damage
                SOURCES     ${Damage_sources}
                HEADERS     ${Damage_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Damage
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Damage_headers}
        DESTINATION include/Damage
        )
