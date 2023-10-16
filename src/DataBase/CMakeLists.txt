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
    FieldUpdatePolicy.hh
    FieldUpdatePolicyInline.hh
    CopyState.hh
    CopyStateInline.hh
    IncrementBoundedState.hh
    IncrementBoundedStateInline.hh
    IncrementState.hh
    IncrementStateInline.hh
    ReplaceBoundedState.hh
    ReplaceBoundedStateInline.hh
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

spheral_add_obj_library(DataBase SPHERAL_OBJ_LIBS)

