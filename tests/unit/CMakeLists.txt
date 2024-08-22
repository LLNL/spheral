add_subdirectory(CXXTests)

if (ENABLE_TIMER)
  set(MPI_TIMER_VAR "True")
else()
  set(MPI_TIMER_VAR "False")
endif()
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Utilities/testTimers.py.in"
  "${SPHERAL_TEST_INSTALL_PREFIX}/tests/unit/Utilities/testTimers.py"
  )
