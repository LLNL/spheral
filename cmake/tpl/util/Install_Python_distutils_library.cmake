# General rule for installing pip modules. Downloads necessary tar/wheel files to the cache dir
# then installs in the future from these files so that offline builds are possible.
#
# NOTE:
# Since we are still on Python 2.7.*, we have to make a few compromises here:
#   1. We lock the exact versions of PIP packages we're using to the older known Python 2 friendly versions.
#   2. We removed the --upgrade option from the PIP install command.
# Hopefully once we move to Python 3 we can relax this stringent version locking and just stay up to date.

# We want to expose python lib directories to the cmake command line. This helps with building Spheral
# with externally installed python modules such as those created by Spack.

macro(Install_Pip_Module lib_name)

  string(REGEX REPLACE "[\=]+[^ ]*" "" lib_name_str ${lib_name})

  set(${lib_name_str}_TARGET_DEPENDS )
  set(${lib_name_str}_download_stamp_file "${CACHE_DIR}/.${lib_name_str}_pip_download.stamp")
  add_custom_command(
    OUTPUT ${${lib_name_str}_download_stamp_file}
    COMMAND echo "-- pip downloading ${lib_name}"
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${BUILDTIME_PYTHONENV_STR} ${PYTHON_EXE} ${PIP_EXE} ${OUT_PROTOCOL_PIP} download --no-binary :all -d ${CACHE_DIR} ${lib_name}
    COMMAND touch ${${lib_name_str}_download_stamp_file}
  )

  set(commands )
  list(APPEND commands COMMAND echo "-- pip install ${lib_name_str}")
  list(APPEND commands COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${BUILDTIME_PYTHONENV_STR} ${PYTHON_EXE} ${PIP_EXE} ${OUT_PROTOCOL_PIP} install ${lib_name} --no-index --find-links ${CACHE_DIR} --target ${${lib_name_str}_DIR})

  set(${lib_name_str}_install_stamp_file "${CACHE_DIR}/.${lib_name_str}_pip_install.stamp")
  add_custom_command(
    OUTPUT ${${lib_name_str}_install_stamp_file}
    ${commands}
    COMMAND touch ${${lib_name_str}_install_stamp_file}
    DEPENDS pip-install ${${lib_name_str}_DEPENDS} ${${lib_name_str}_download_stamp_file}
  )

  set(${lib_name_str}_TARGET_DEPENDS ${${lib_name_str}_download_stamp_file} ${${lib_name_str}_install_stamp_file})
  
  add_custom_target(
    ${lib_name_str}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${BUILDTIME_PYTHONENV_STR} ${PYTHON_EXE} -c \"import ${${lib_name_str}_IMPORT}\" 
    DEPENDS pip-install ${${lib_name_str}_TARGET_DEPENDS} 
    )

endmacro()

