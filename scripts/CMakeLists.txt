configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/spheral-setup-venv.in"
  "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/spheral-env.in"
  "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
  )

install(FILES 
    "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
    "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts"
  )

install(CODE "execute_process( \
  COMMAND bash ./scripts/spheral-setup-venv.sh \ 
  )")
