#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_spheral_tests
   )


#------------------------------------------------------------------------------
# Add gtest C++ tests
#------------------------------------------------------------------------------

get_property(SPHERAL_CXX_LIBS GLOBAL PROPERTY SPHERAL_CXX_LIBS)

include_directories(${SPHERAL_ROOT_DIR}/src)

set(TESTS)

foreach(test ${gtest_spheral_tests})
  get_filename_component( test_name ${test} NAME_WE )

  list(APPEND TESTS '${test_name}.py')

  blt_add_executable(NAME ${test_name}
                     SOURCES ${test}
                     OUTPUT_DIR ${PROJECT_BINARY_DIR}/tests/unit/CXXTests
                     DEPENDS_ON gtest -Wl,--start-group ${SPHERAL_CXX_LIBS} ${spheral_blt_depends} -Wl,--end-group ${CMAKE_DL_LIBS}
                     )

  blt_add_test(NAME ${test_name}
               COMMAND ${test_name}
               )

  if (NOT ENABLE_CXXONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pyRunCXXTest.in"
                   "${SPHERAL_TEST_INSTALL_PREFIX}/unit/CXXTests/${test_name}.py"
                   )
  endif()
  
  install(TARGETS ${test_name}
          DESTINATION ${SPHERAL_TEST_INSTALL_PREFIX}/tests/unit/CXXTests
          )
  
  set_target_properties(${test_name} PROPERTIES
                        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${qhull_DIR}/lib;${conduit_DIR}/lib;${axom_DIR}/lib;${boost_DIR}/lib"
                        )
endforeach()

string(REPLACE ";" ", " TEST_LIST "${TESTS}")

if (NOT ENABLE_CXXONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/runCXXTests.in"
                 "${SPHERAL_TEST_INSTALL_PREFIX}/unit/CXXTests/runCXXTests.ats"
                 )
endif()
