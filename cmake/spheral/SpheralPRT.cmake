#----------------------------------------------------------------------------------------
#                                   Spheral_Python_Runtime_Env
#----------------------------------------------------------------------------------------

function(Spheral_Python_Env target_name)

  # Define our arguments
  set(options )
  set(oneValueArgs PREFIX)
  set(multiValueArgs REQUIREMENTS)
  cmake_parse_arguments(${target_name} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(REQUIREMENTS_ARGS)
  foreach(_req ${${target_name}_REQUIREMENTS})
    list(APPEND REQUIREMENTS_ARGS -r)
    list(APPEND REQUIREMENTS_ARGS ${SPHERAL_ROOT_DIR}/scripts/${_req})
  endforeach()

  add_custom_target(${target_name} ALL
    COMMAND ${Python3_EXECUTABLE} -m venv ${${target_name}_PREFIX}/.venv;
    COMMAND . ${${target_name}_PREFIX}/.venv/bin/activate &&
            python -m pip install --upgrade pip &&
            python -m pip install ${REQUIREMENTS_ARGS}
    VERBATIM
    DEPENDS Python3::Python
  )
  set_property(TARGET ${target_name} PROPERTY EXECUTABLE python)
  set_property(TARGET ${target_name} PROPERTY ACTIVATE_VENV . ${${target_name}_PREFIX}/.venv/bin/activate)
endfunction()

