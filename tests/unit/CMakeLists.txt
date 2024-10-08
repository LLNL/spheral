add_subdirectory(CXXTests)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/testTimers.py.in"
  "${SPHERAL_TEST_INSTALL_PREFIX}/tests/unit/Utilities/testTimers.py"
  )
