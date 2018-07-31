set(Damage_sources2d
    DamageModelInst.cc.2.cc
    JohnsonCookDamageInst.cc.2.cc
    JohnsonCookDamagePolicyInst.cc.2.cc
    JohnsonCookFailureStrainPolicyInst.cc.2.cc
    LongitudinalSoundSpeedPolicyInst.cc.2.cc
    StrainPolicyInst.cc.2.cc
    TensorDamageModelInst.cc.2.cc
    TensorDamagePolicyInst.cc.2.cc
    TensorStrainPolicyInst.cc.2.cc
    YoungsModulusPolicyInst.cc.2.cc
    computeFragmentFieldInst.cc.2.cc
    weibullFlawDistributionInst.cc.2.cc
    )

set(Damage_sources3d
    DamageModelInst.cc.3.cc
    JohnsonCookDamageInst.cc.3.cc
    JohnsonCookDamagePolicyInst.cc.3.cc
    JohnsonCookFailureStrainPolicyInst.cc.3.cc
    LongitudinalSoundSpeedPolicyInst.cc.3.cc
    StrainPolicyInst.cc.3.cc
    TensorDamageModelInst.cc.3.cc
    TensorDamagePolicyInst.cc.3.cc
    TensorStrainPolicyInst.cc.3.cc
    YoungsModulusPolicyInst.cc.3.cc
    computeFragmentFieldInst.cc.3.cc
    weibullFlawDistributionInst.cc.3.cc
    )

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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Damage
                SOURCES     ${Damage_sources2d} ${Damage_sources3d}
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
