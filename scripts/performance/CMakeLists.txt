
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/performance/performance.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/performance/performance.py"
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/performance/performance.py"
  "${CMAKE_CURRENT_BINARY_DIR}/performance/test_list.json"
  "${CMAKE_CURRENT_BINARY_DIR}/performance/compare_times.py"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/tests/performance"
)
