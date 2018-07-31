set(FieldOperations_sources2d
    binFieldList2LatticeInst.cc.2.cc
    divergenceInst.cc.2.cc
    divergenceMashInst.cc.2.cc
    divergencePairWiseInst.cc.2.cc
    gradDivVectorFieldListInst.cc.2.cc
    gradDivVectorFieldListGoldenInst.cc.2.cc
    gradDivVectorFieldListGolden2Inst.cc.2.cc
    gradDivVectorFieldListMashInst.cc.2.cc
    gradDivVectorFieldListPairWiseInst.cc.2.cc
    gradDivVectorFieldListSimpleInst.cc.2.cc
    gradientInst.cc.2.cc
    gradientMashInst.cc.2.cc
    gradientPairWiseInst.cc.2.cc
    sampleFieldsMashInst.cc.2.cc
    sampleMultipleFields2LatticeInst.cc.2.cc
    sampleMultipleFields2LatticeMashInst.cc.2.cc
    sampleMultipleFieldsMashInst.cc.2.cc
    smoothFieldsInst.cc.2.cc
    smoothFieldsMashInst.cc.2.cc
    smoothFieldsMash2Inst.cc.2.cc
    splatFieldsMashInst.cc.2.cc
    splatMultipleFieldsMashInst.cc.2.cc
    )

set(FieldOperations_sources3d
    binFieldList2LatticeInst.cc.3.cc
    divergenceInst.cc.3.cc
    divergenceMashInst.cc.3.cc
    divergencePairWiseInst.cc.3.cc
    gradDivVectorFieldListInst.cc.3.cc
    gradDivVectorFieldListGoldenInst.cc.3.cc
    gradDivVectorFieldListGolden2Inst.cc.3.cc
    gradDivVectorFieldListMashInst.cc.3.cc
    gradDivVectorFieldListPairWiseInst.cc.3.cc
    gradDivVectorFieldListSimpleInst.cc.3.cc
    gradientInst.cc.3.cc
    gradientMashInst.cc.3.cc
    gradientPairWiseInst.cc.3.cc
    sampleFieldsMashInst.cc.3.cc
    sampleMultipleFields2LatticeInst.cc.3.cc
    sampleMultipleFields2LatticeMashInst.cc.3.cc
    sampleMultipleFieldsMashInst.cc.3.cc
    smoothFieldsInst.cc.3.cc
    smoothFieldsMashInst.cc.3.cc
    smoothFieldsMash2Inst.cc.3.cc
    splatFieldsMashInst.cc.3.cc
    splatMultipleFieldsMashInst.cc.3.cc
    )

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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_FieldOperations
                SOURCES     ${FieldOperations_sources2d} ${FieldOperations_sources3d}
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
