set(Field_sources2d
    AllNodeIteratorInst.cc.2.cc
    CoarseNodeIteratorInst.cc.2.cc
    GhostNodeIteratorInst.cc.2.cc
    InternalNodeIteratorInst.cc.2.cc
    MasterNodeIteratorInst.cc.2.cc
    NodeIteratorBaseInst.cc.2.cc
    RefineNodeIteratorInst.cc.2.cc
    )

set(Field_sources3d
    AllNodeIteratorInst.cc.3.cc
    CoarseNodeIteratorInst.cc.3.cc
    GhostNodeIteratorInst.cc.3.cc
    InternalNodeIteratorInst.cc.3.cc
    MasterNodeIteratorInst.cc.3.cc
    NodeIteratorBaseInst.cc.3.cc
    RefineNodeIteratorInst.cc.3.cc
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

INCLUDE_DIRECTORIES(../)

set(spheral_depends BOOST Eigen qhull silo polytope mpi)

blt_add_library(NAME        Spheral_Field
                SOURCES     ${Field_sources2d} ${Field_sources3d}
                HEADERS     ${Field_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_Field
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Field_headers}
        DESTINATION include/Field
        )
