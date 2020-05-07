if (NOT TPL_PARALLEL_BUILD)
  set (TPL_PARALLEL_BUILD "1")
endif()

if (TPL_VERBOSE)
  set(OUT_PROTOCOL_EP 0)
else()
  set(OUT_PROTOCOL_EP 1)
  set(OUT_PROTOCOL_PIP "-q")
endif()

if (NOT SPHERAL_INSTALL_DIR)
  get_filename_component(DEFAULT_TPL_LOCATION ${CMAKE_BINARY_DIR}/Spheral-tpl/ ABSOLUTE)
else()
  set(DEFAULT_TPL_LOCATION ${SPHERAL_INSTALL_DIR})
endif()
message("Default TPL location : ${DEFAULT_TPL_LOCATION}\n")

#----------------------------------------------------------------------------------------
#                                   Spheral_Handle_TPL
#----------------------------------------------------------------------------------------

# -------------------------------------------
# VARIBALES THAT NEED TO BE PREVIOSLY DEFINED
# -------------------------------------------
# BUILD_TPL             : REQUIRED : Set at configure time defines weather to search or build tpl.
# <lib_name>_DIR        : OPTIONAL : The installation location of the tpl, if installed or not.
#                                    If not set and BUILD_TPL=On <lib_name> will be installed to a 
#                                    default loc. 

# ----------------------
# INPUT-OUTPUT VARIBALES
# ----------------------
# <lib_name>   : REQUIRED : The name of the target tpl
# <dep_list>   : REQUIRED : list that appends new target <lib_name> to itself.

# -----------------------
# OUTPUT VARIABLES TO USE - Made available implicitly after funciton call.
# -----------------------
# <lib_name>_libs : list of full paths to tpl lib files to be linked to.

function(Spheral_Handle_TPL lib_name dep_list)
                                                       
  if(NOT BUILD_TPL OR NOT ${lib_name}_BUILD)

    if (NOT ${lib_name}_DIR)
      set(${lib_name}_DIR ${DEFAULT_TPL_LOCATION}/${lib_name})
      message("${lib_name}_DIR not set.")
      message("Setting ${lib_name} search to default location : ${${lib_name}_DIR}")
    else()
      message("${lib_name}_DIR set.")
      message("Searching ${lib_name} for : ${${lib_name}_DIR}")
    endif()

    if(NOT BUILD_TPL)
      set(${lib_name}_BUILD Off)
    endif()

  else()

    message("Generating build process for ${lib_name}.")

    if (NOT ${lib_name}_DIR)
      set(${lib_name}_DIR ${DEFAULT_TPL_LOCATION}/${lib_name})
      message("${lib_name}_DIR not set. Installing ${lib_name} to default location : ${${lib_name}_DIR}")
    else()
      message("${lib_name}_DIR set. Installing ${lib_name} to : ${${lib_name}_DIR}")
    endif()
  endif()

  set(${lib_name}_ADD_BLT_TARGET ON)
  include(${TPL_CMAKE_DIR}/${lib_name}.cmake)

  list(APPEND ${lib_name}_INCLUDES $<BUILD_INTERFACE:${${lib_name}_DIR}/include>)

  # Generate full path to lib file for output list.
  set(${lib_name}_LIBRARIES )
  foreach(lib ${${lib_name}_libs})
    set(lib_abs "${${lib_name}_DIR}/lib/${lib}")
    list(APPEND ${lib_name}_LIBRARIES $<BUILD_INTERFACE:${lib_abs}>)

    # Check all necessary files exist during build time when not installing TPL
    if (NOT BUILD_TPL OR NOT ${lib_name}_BUILD)
      if (NOT EXISTS ${lib_abs})
        message(FATAL_ERROR "Cannot find ${lib} in ${${lib_name}_DIR} for TPL ${lib_name}.")
      else()
        message("Found: ${lib_abs}")
      endif()
    endif()
  endforeach()

  blt_register_library(
    NAME blt_${lib_name}
    INCLUDES ${${lib_name}_INCLUDES}
    LIBRARIES ${${lib_name}_LIBRARIES}
    )

  if (${lib_name}_ADD_BLT_TARGET)
    list(APPEND spheral_blt_depends blt_${lib_name})
    if (${lib_name}_BUILD)
      list(APPEND ${dep_list} ${lib_name})
    endif()
  endif()

  set(${dep_list} ${${dep_list}} PARENT_SCOPE)
  set(spheral_blt_depends ${spheral_blt_depends} PARENT_SCOPE)

  message("")

endfunction()
