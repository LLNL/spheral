if (NOT ENABLE_CXXONLY)

  set(SPHERAL_ATS_BUILD_CONFIG_ARGS )
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"level<100\"'")
  endif()
  if (NOT ENABLE_MPI)
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"np<2\"'")
  endif()
  if (ENABLE_CUDA)
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--addOp --smpi_off")
  endif()

  string(REPLACE ";" " " SPHERAL_ATS_BUILD_CONFIG_ARGS_STRING "${SPHERAL_ATS_BUILD_CONFIG_ARGS}")
  
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/spheral-setup-venv.in"
    "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
    )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/spheral-env.in"
    "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
    )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/atstest.in"
    "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
    )

  install(FILES 
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts"
    )

  install(CODE "execute_process( \
    COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh \ 
    )")
endif()
