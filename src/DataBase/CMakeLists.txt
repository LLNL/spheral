set(DataBase_sources2d
    CompositeFieldListPolicyInst.cc.2.cc
    CopyFieldListInst.cc.2.cc
    CopyStateInst.cc.2.cc
    DataBaseInst.cc.2.cc
    IncrementBoundedFieldListInst.cc.2.cc
    IncrementBoundedStateInst.cc.2.cc
    IncrementFieldListInst.cc.2.cc
    IncrementStateInst.cc.2.cc
    ReplaceBoundedFieldListInst.cc.2.cc
    ReplaceBoundedStateInst.cc.2.cc
    ReplaceFieldListInst.cc.2.cc
    ReplaceStateInst.cc.2.cc
    StateInst.cc.2.cc
    StateBaseInst.cc.2.cc
    StateDerivativesInst.cc.2.cc
    )

set(DataBase_sources3d
    CompositeFieldListPolicyInst.cc.3.cc
    CopyFieldListInst.cc.3.cc
    CopyStateInst.cc.3.cc
    DataBaseInst.cc.3.cc
    IncrementBoundedFieldListInst.cc.3.cc
    IncrementBoundedStateInst.cc.3.cc
    IncrementFieldListInst.cc.3.cc
    IncrementStateInst.cc.3.cc
    ReplaceBoundedFieldListInst.cc.3.cc
    ReplaceBoundedStateInst.cc.3.cc
    ReplaceFieldListInst.cc.3.cc
    ReplaceStateInst.cc.3.cc
    StateInst.cc.3.cc
    StateBaseInst.cc.3.cc
    StateDerivativesInst.cc.3.cc
    )

set(DataBase_headers
    CompositeFieldListPolicy.hh
    CopyFieldList.hh
    CopyState.hh
    DataBase.hh
    DataBaseInline.hh
    FieldListUpdatePolicyBase.hh
    FieldListUpdatePolicyBaseInline.hh
    FieldUpdatePolicyBase.hh
    FieldUpdatePolicyBaseInline.hh
    IncrementBoundedFieldList.hh
    IncrementBoundedFieldListInline.hh
    IncrementBoundedState.hh
    IncrementBoundedStateInline.hh
    IncrementFieldList.hh
    IncrementState.hh
    ReplaceBoundedFieldList.hh
    ReplaceBoundedFieldListInline.hh
    ReplaceBoundedState.hh
    ReplaceBoundedStateInline.hh
    ReplaceFieldList.hh
    ReplaceState.hh
    State.hh
    StateBase.hh
    StateBaseInline.hh
    StateDerivatives.hh
    StateDerivativesInline.hh
    StateInline.hh
    UpdatePolicyBase.hh
    UpdatePolicyBaseInline.hh
    )

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_DataBase
                SOURCES     ${DataBase_sources2d} ${DataBase_sources3d}
                HEADERS     ${DataBase_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_DataBase
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${DataBase_headers}
        DESTINATION include/DataBase
        )
