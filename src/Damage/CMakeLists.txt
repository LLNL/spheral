include_directories(.)

set(Damage_sources )

# Add explicit instantiations
list(APPEND Damage_sources
     DamageModelInst.cc
     TensorDamageModelInst.cc
     ProbabilisticDamageModelInst.cc
     ProbabilisticDamagePolicyInst.cc
     IvanoviSALEDamageModelInst.cc
     IvanoviSALEDamagePolicyInst.cc
     StrainPolicyInst.cc
     TensorStrainPolicyInst.cc
     TensorDamagePolicyInst.cc
     YoungsModulusPolicyInst.cc
     LongitudinalSoundSpeedPolicyInst.cc
     weibullFlawDistributionBenzAsphaugInst.cc
     weibullFlawDistributionOwenInst.cc
     computeFragmentFieldInst.cc
     JohnsonCookDamageInst.cc
     JohnsonCookFailureStrainPolicyInst.cc
     JohnsonCookDamagePolicyInst.cc
     PairMaxDamageNodeCouplingInst.cc
     ThreePointDamagedNodeCouplingInst.cc
     DamageGradientNodeCouplingInst.cc
     )

set(Damage_headers
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
    ScalarDamageModel.hh
    ScalarDamagePolicy.hh
    StrainPolicy.hh
    TensorDamageModel.hh
    TensorDamageModelInline.hh
    ProbabilisticDamageModel.hh
    ProbabilisticDamageModelInline.hh
    IvanoviSALEDamageModel.hh
    IvanoviSALEDamageModelInline.hh
    TensorDamagePolicy.hh
    ProbabilisticDamagePolicy.hh
    IvanoviSALEDamagePolicy.hh
    TensorStrainPolicy.hh
    YoungsModulusPolicy.hh
    computeFragmentField.hh
    oneMinusDamage.hh
    weibullFlawDistributionBenzAsphaug.hh
    weibullFlawDistributionOwen.hh
    DamageGradientNodeCoupling.hh
    PairMaxDamageNodeCoupling.hh
    MaxDamageNodeCoupling.hh
    ThreePointDamagedNodeCoupling.hh
    )

spheral_install_python_files(
  identifyFragments.py
  WeibullDamage.py
  GradyKippTensorDamage.py
  JohnsonCookDamageFactories.py
  ShadowProbabilisticDamageModel.py
  ShadowIvanoviSALEDamageModel.py
  )

spheral_add_obj_library(Damage SPHERAL_OBJ_LIBS)

