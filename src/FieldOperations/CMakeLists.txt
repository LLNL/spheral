include_directories(.)

set(FieldOperations_sources )

# Add explicit instantiations
list(APPEND FieldOperations_sources
    smoothFieldsInst.cc
    gradientInst.cc
    divergenceInst.cc
    smoothFieldsMashInst.cc
    smoothFieldsMash2Inst.cc
    sampleFieldsMashInst.cc
    splatFieldsMashInst.cc
    gradientMashInst.cc
    divergenceMashInst.cc
    gradientPairWiseInst.cc
    divergencePairWiseInst.cc
    sampleMultipleFieldsMashInst.cc
    splatMultipleFieldsMashInst.cc
    gradDivVectorFieldListInst.cc
    gradDivVectorFieldListSimpleInst.cc
    gradDivVectorFieldListGoldenInst.cc
    gradDivVectorFieldListGolden2Inst.cc
    gradDivVectorFieldListMashInst.cc
    gradDivVectorFieldListPairWiseInst.cc
    sampleMultipleFields2LatticeInst.cc
    sampleMultipleFields2LatticeMashInst.cc
    binFieldList2LatticeInst.cc
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

spheral_add_obj_library(FieldOperations SPHERAL_OBJ_LIBS)

