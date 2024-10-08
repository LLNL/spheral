
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/performance.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/performance.py"
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/performance.py"
  DESTINATION "${SPHERAL_TEST_INSTALL_PREFIX}"
)
