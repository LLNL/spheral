set(CXXTests_inst
  #testNodeIterators
  )

set(CXXTests_sources
  #test_r3d_utils.cc
  #testNodeIterators.cc
  )

instantiate(CXXTests_inst CXXTests_sources)

set(CXXTests_headers
  #test_r3d_utils.hh
  testFieldListIndexing.hh
  testNewtonRaphson.hh
  #testNodeIterators.hh
  )

set(CXXTests_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/CXXTests
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_CXXTests
                HEADERS     ${CXXTests_headers}
                SOURCES     ${CXXTests_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${CXXTests_includes}
                )

install(TARGETS             Spheral_CXXTests
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${CXXTests_headers}
        DESTINATION include/CXXTests
        )

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_CXXTests)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
