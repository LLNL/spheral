set(Damage_inst
  computeFragmentField
  DamageModel
  JohnsonCookDamage
  JohnsonCookDamagePolicy
  JohnsonCookFailureStrainPolicy
  LongitudinalSoundSpeedPolicy
  StrainPolicy
  TensorDamageModel
  TensorDamagePolicy
  TensorStrainPolicy
  weibullFlawDistribution
  YoungsModulusPolicy
  )

set(Damage_sources
  #DamagedSoundSpeedPolicy.cc
  #GradyKippScalarDamage.cc
  #ScalarDamageModel.cc
  #ScalarDamagePolicy.cc
  )
instantiate(Damage_inst Damage_sources)

set(Damage_headers
  computeFragmentField.hh
  #DamagedSoundSpeedPolicy.hh
  DamageGradientPolicy.hh
  DamageModel.hh
  DamageModelInline.hh
  EffectiveTensorDamagePolicy.hh
  #GradyKippScalarDamage.hh
  JohnsonCookDamage.hh
  JohnsonCookDamageInline.hh
  JohnsonCookDamagePolicy.hh
  JohnsonCookFailureStrainPolicy.hh
  LongitudinalSoundSpeedPolicy.hh
  oneMinusDamage.hh
  #ScalarDamageModel.hh
  #ScalarDamagePolicy.hh
  StrainPolicy.hh
  TensorDamageModel.hh
  TensorDamageModelInline.hh
  TensorDamagePolicy.hh
  TensorStrainPolicy.hh
  WeibullDamage.py
  weibullFlawDistribution.hh
  YoungsModulusPolicy.hh
  )

set(Damage_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Damage
  ${SPHERAL_SOURCE_DIR}/Strength
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_Damage
                HEADERS     ${Damage_headers}
                SOURCES     ${Damage_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Damage_includes}
                )

install(TARGETS             Spheral_Damage
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Damage_headers}
        DESTINATION include/Damage
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Damage)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
