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

set(Field_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Field
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_Field
                HEADERS     ${Field_headers}
                SOURCES     ${Field_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Field_includes}
                )

install(TARGETS             Spheral_Field
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Field_headers}
        DESTINATION include/Field
        )
