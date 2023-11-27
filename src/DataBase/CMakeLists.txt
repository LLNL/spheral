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
    DataBase.hh
    DataBaseInline.hh
    UpdatePolicyBase.hh
    UpdatePolicyBaseInline.hh
    FieldUpdatePolicy.hh
    FieldUpdatePolicyInline.hh
    CopyState.hh
    CopyStateInline.hh
    IncrementState.hh
    IncrementStateInline.hh
    IncrementBoundedState.hh
    IncrementBoundedStateInline.hh
    ReplaceState.hh
    ReplaceStateInline.hh
    ReplaceBoundedState.hh
    ReplaceBoundedStateInline.hh
    State.hh
    StateBase.hh
    StateBaseInline.hh
    StateDerivatives.hh
    StateDerivativesInline.hh
    StateInline.hh
    )

spheral_add_obj_library(DataBase SPHERAL_OBJ_LIBS)

