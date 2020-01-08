include_directories(.)
set(Field_inst
    AllNodeIterator
    CoarseNodeIterator
    GhostNodeIterator
    InternalNodeIterator
    MasterNodeIterator
    NodeIteratorBase
    RefineNodeIterator
   )


set(Field_sources )

instantiate(Field_inst Field_sources)

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

spheral_add_cxx_library(Field)

install(TARGETS              Spheral_Field
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Field_headers}
        DESTINATION include/Field
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Field)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
