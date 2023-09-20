include_directories(.)
set(DataBase_inst
    DataBase
    State
    StateBase
    StateDerivatives
   )


set(DataBase_sources )

instantiate(DataBase_inst DataBase_sources)

set(DataBase_headers
    CompositeFieldListPolicy.hh
    CompositeFieldListPolicyInline.hh
    CopyFieldList.hh
    CopyFieldListInline.hh
    CopyState.hh
    CopyStateInline.hh
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
    IncrementFieldListInline.hh
    IncrementState.hh
    IncrementStateInline.hh
    ReplaceBoundedFieldList.hh
    ReplaceBoundedFieldListInline.hh
    ReplaceBoundedState.hh
    ReplaceBoundedStateInline.hh
    ReplaceFieldList.hh
    ReplaceFieldListInline.hh
    ReplaceState.hh
    ReplaceStateInline.hh
    State.hh
    StateBase.hh
    StateBaseInline.hh
    StateDerivatives.hh
    StateDerivativesInline.hh
    StateInline.hh
    UpdatePolicyBase.hh
    UpdatePolicyBaseInline.hh
    )

spheral_add_obj_library(DataBase)

