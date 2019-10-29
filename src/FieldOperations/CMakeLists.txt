set(FieldOperations_inst
  binFieldList2Lattice
  divergence
  divergenceMash
  divergencePairWise
  gradDivVectorFieldListGolden2
  gradDivVectorFieldListGolden
  gradDivVectorFieldList
  gradDivVectorFieldListMash
  gradDivVectorFieldListPairWise
  gradDivVectorFieldListSimple
  gradient
  gradientMash
  gradientPairWise
  sampleFieldsMash
  sampleMultipleFields2Lattice
  sampleMultipleFields2LatticeMash
  sampleMultipleFieldsMash
  smoothFields
  smoothFieldsMash2
  smoothFieldsMash
  splatFieldsMash
  splatMultipleFieldsMash
  )

set(FieldOperations_sources )
instantiate(FieldOperations_inst FieldOperations_sources)

set(FieldOperations_headers
  binFieldList2Lattice.hh
  FieldListFunctions.hh
  FieldListFunctionsMash.hh
  FieldListSecondDerivatives.hh
  incrementPairWiseDelta.hh
  monotonicallyLimitedGradient.hh
  PairWiseFieldListFunctions.hh
  sampleMultipleFields2Lattice.hh
  )

set(FieldOperations_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/FieldOperations
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_FieldOperations
                HEADERS     ${FieldOperations_headers}
                SOURCES     ${FieldOperations_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${FieldOperations_includes}
                )

install(TARGETS             Spheral_FieldOperations
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${FieldOperations_headers}
        DESTINATION include/FieldOperations
        )
