#-----------------------------------------------------------------------------------
# spheral_install_python_files
#     - Copies Python files to the Spheral install path, and byte compiles them
#
# pyfiles : a list of python files to install
#
# Note, if ENABLE_CXXONLY is set, this function does nothing
#-----------------------------------------------------------------------------------

function(spheral_install_python_files pyfiles)

  if (NOT ENABLE_CXXONLY)
    install(FILES ${pyfiles}
            DESTINATION Spheral)
    install(CODE "execute_process( \
            COMMAND ${SPHERAL_INSTALL_DIR}/python/bin/python -m compileall Spheral \
            WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX})")
  endif()

endfunction()
