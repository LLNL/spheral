if (NOT ENABLE_CXXONLY)

  set(SPHERAL_ATS_BUILD_CONFIG_ARGS )

  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"level<100\"'")
  endif()

  if (NOT ENABLE_MPI)
    list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"np<2\"'")
  endif()

  if (NOT SPHERAL_ENABLE_FSISPH)
     list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"not fsisph\"'")
  endif()

  if (NOT SPHERAL_ENABLE_GSPH)
     list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"not gsph\"'")
  endif()

  if (NOT SPHERAL_ENABLE_SVPH)
     list(APPEND SPHERAL_ATS_BUILD_CONFIG_ARGS "--filter='\"not svph\"'")
  endif()

  if ($ENV{SYS_TYPE} MATCHES ".*blueos.*")
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

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/lcatstest.in"
    "${CMAKE_CURRENT_BINARY_DIR}/lcatstest.sh"
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/performance/performance.py.in"
    "${CMAKE_CURRENT_BINARY_DIR}/performance/performance.py"
  )

  install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/spheral-env.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/atstest.sh"
      "${CMAKE_CURRENT_BINARY_DIR}/lcatstest.sh"
      "${CMAKE_CURRENT_SOURCE_DIR}/lc/lcats"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/scripts"
  )

  install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/performance/performance.py"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/tests"
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/runtime-requirements.txt.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/runtime-requirements.txt"
  )

  # Install Runtime python dependencies into the Install directory.
  Spheral_Python_Env(python_runtime_env 
    REQUIREMENTS runtime-requirements.txt 
    PREFIX ${CMAKE_INSTALL_PREFIX}
  )

  install(FILES ${polytope_DIR}/${SPHERAL_SITE_PACKAGES_PATH}/polytope/polytope.so
    DESTINATION ${CMAKE_INSTALL_PREFIX}/.venv/${SPHERAL_SITE_PACKAGES_PATH}/polytope/
  )

  install(CODE "execute_process( \
    COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/spheral-setup-venv.sh)"
  )
endif()
