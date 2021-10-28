#------------------------------------------------------------------------------
# Specify list of tests
#------------------------------------------------------------------------------

set(gtest_spheral_tests
    testSidreStorage.cc
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
                     OUTPUT_DIR ${PROJECT_BINARY_DIR}/tests/Unit/CXXTests
                     DEPENDS_ON gtest -Wl,--start-group ${SPHERAL_CXX_LIBS} ${spheral_blt_depends} -Wl,--end-group ${CMAKE_DL_LIBS}
                     )

  blt_add_test(NAME ${test_name}
               COMMAND ${test_name}
               )

  if (NOT ENABLE_CXXONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/pyRunCXXTest.in"
                   "${CMAKE_INSTALL_PREFIX}/tests/Unit/CXXTests/${test_name}.py"
                   )
  endif()
  
  install(TARGETS ${test_name}
          DESTINATION ${CMAKE_INSTALL_PREFIX}/tests/Unit/CXXTests
          )
  
  set_target_properties(${test_name} PROPERTIES
                        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;"
                        )
endforeach()

string(REPLACE ";" ", " TEST_LIST "${TESTS}")

if (NOT ENABLE_CXXONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/runCXXTests.in"
                 "${CMAKE_INSTALL_PREFIX}/tests/Unit/CXXTests/runCXXTests.ats"
                 )
endif()
