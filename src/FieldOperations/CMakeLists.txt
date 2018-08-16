set(FieldOperations_inst
    binFieldList2Lattice
    divergence
    divergenceMash
    divergencePairWise
    gradDivVectorFieldList
    gradDivVectorFieldListGolden
    gradDivVectorFieldListGolden2
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
    smoothFieldsMash
    smoothFieldsMash2
    splatFieldsMash
    splatMultipleFieldsMash
   )


set(FieldOperations_sources )

instantiate(FieldOperations_inst FieldOperations_sources)

set(FieldOperations_headers
    FieldListFunctions.hh
    FieldListFunctionsMash.hh
    FieldListSecondDerivatives.hh
    PairWiseFieldListFunctions.hh
    binFieldList2Lattice.hh
    incrementPairWiseDelta.hh
    monotonicallyLimitedGradient.hh
    sampleMultipleFields2Lattice.hh
    )

INCLUDE_DIRECTORIES(. ../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_FieldOperations
                SOURCES     ${FieldOperations_sources}
                HEADERS     ${FieldOperations_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_FieldOperations
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${FieldOperations_headers}
        DESTINATION include/FieldOperations
        )
