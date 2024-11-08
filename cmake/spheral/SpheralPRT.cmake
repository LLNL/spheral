#----------------------------------------------------------------------------------------
#                                   Spheral_Python_Runtime_Env
#----------------------------------------------------------------------------------------

function(Spheral_Python_Env target_name requirements_file prefix)
  add_custom_target(${target_name} ALL
    COMMAND ${Python3_EXECUTABLE} -m venv ${prefix}/.venv;
    COMMAND . ${prefix}/.venv/bin/activate &&
            python -m pip install --upgrade pip &&
            python -m pip install -r ${SPHERAL_ROOT_DIR}/scripts/${requirements_file}
    VERBATIM
    DEPENDS Python3::Python
  )
  set_property(TARGET ${target_name} PROPERTY EXECUTABLE python)
  set_property(TARGET ${target_name} PROPERTY ACTIVATE_VENV . ${prefix}/.venv/bin/activate)
endfunction()

