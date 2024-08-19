include_directories(.)

set(Field_sources )

# Add explicit instantiations
list(APPEND Field_sources
     AllNodeIteratorInst.cc
     CoarseNodeIteratorInst.cc
     GhostNodeIteratorInst.cc
     InternalNodeIteratorInst.cc
     MasterNodeIteratorInst.cc
     NodeIteratorBaseInst.cc
     RefineNodeIteratorInst.cc
     )

set(Field_headers
    AggregateField.hh
    AllNodeIterator.hh
    AllNodeIteratorInline.hh
    CoarseNodeIterator.hh
    CoarseNodeIteratorInline.hh
    Field.hh
    FieldBase.hh
    FieldBaseInline.hh
    FieldInline.hh
    FieldList.hh
    FieldListBase.hh
    FieldListBaseInline.hh
    FieldListInline.hh
    FieldListSet.hh
    FieldListSetInline.hh
    FieldListThreadWrapper.hh
    GhostNodeIterator.hh
    GhostNodeIteratorInline.hh
    InternalNodeIterator.hh
    InternalNodeIteratorInline.hh
    MasterNodeIterator.hh
    MasterNodeIteratorInline.hh
    NodeIteratorBase.hh
    NodeIteratorBaseInline.hh
    NodeIterators.hh
    RefineNodeIterator.hh
    RefineNodeIteratorInline.hh
    uvm_allocator.hh
    )

spheral_add_obj_library(Field SPHERAL_OBJ_LIBS)

