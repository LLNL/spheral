#-------------------------------------------------------------------------------'
# BLT targets don't play well with standard CMake targets, so we extract the 
# necessary bits ourselves.
#-------------------------------------------------------------------------------'
function(spheral_extract_blt_props package_name)
  set(${package_name}_includes)
  set(${package_name}_libs)
  string(TOUPPER ${package_name} _target_upper)
  set(_target_prefix "_BLT_${_target_upper}_")
  get_cmake_property(_variable_names VARIABLES)
  foreach (prop ${_variable_names})
    if(prop MATCHES "^${_target_prefix}_INCLUDES")
      list(APPEND ${package_name}_includes ${${prop}})
    elseif(prop MATCHES "${_target_prefix}_LIBRARIES")
      list(APPEND ${package_name}_libs ${${prop}})
    endif()
  endforeach()
  set(${package_name}_includes ${package_name}_includes PARENT_SCOPE)
  set(${package_name}_libs ${package_name}_libs PARENT_SCOPE)
endfunction()
