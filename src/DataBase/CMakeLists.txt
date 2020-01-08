include_directories(.)
set(DataBase_inst
    CompositeFieldListPolicy
    CopyFieldList
    CopyState
    DataBase
    IncrementBoundedFieldList
    IncrementBoundedState
    IncrementFieldList
    IncrementState
    ReplaceBoundedFieldList
    ReplaceBoundedState
    ReplaceFieldList
    ReplaceState
    State
    StateBase
    StateDerivatives
   )


set(DataBase_sources )

instantiate(DataBase_inst DataBase_sources)

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

spheral_add_cxx_library(DataBase)

install(TARGETS              Spheral_DataBase
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${DataBase_headers}
        DESTINATION include/DataBase
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_DataBase)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
