#----------------------------------------------------------------------------------------
#                                   spheral_add_obj_library
#----------------------------------------------------------------------------------------
# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# spheral_blt_depends    : REQUIRED : List of external dependencies
# SPHERAL_CXX_DEPENDS    : REQUIRED : List of compiler dependencies
# <package_name>_headers : OPTIONAL : List of necessary headers to include
# <package_name>_sources : OPTIONAL : List of necessary source files to include
# SPHERAL_SUBMOD_DEPENDS : REQUIRED : List of submodule dependencies
# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# <package_name>   : REQUIRED : Desired package name
# any extra args   : OPTIONAL : Can provide the NAME of a global list
#                               containing additional dependencies
# -----------------------,
# OUTPUT VARIABLES TO USE - Made available implicitly after function call
# -----------------------
# <Spheral_package_name> : Target for a given spheral package
#----------------------------------------------------------------------------------------
function(spheral_add_obj_library package_name)
  # Assumes global variable spheral_blt_depends exists and is filled with external dependencies
  get_property(spheral_blt_depends GLOBAL PROPERTY spheral_blt_depends)
  # Assumes global variable spheral_cxx_depends exists and is filled with compiler dependencies
  get_property(SPHERAL_CXX_DEPENDS GLOBAL PROPERTY SPHERAL_CXX_DEPENDS)
  # For including files in submodules, currently unused
  get_property(SPHERAL_SUBMOD_INCLUDES GLOBAL PROPERTY SPHERAL_SUBMOD_INCLUDES)

  # Only the package_name is expected; any extra arguments are stored in ${ARGN}.
  # This allows the NAME of a global list (not the list itself)
  # with additional dependencies to be passed to the function.

  # First, check if additional arguments are provided
  # CMake requires ${ARGN} is assigned before using with list
  set(extra_args ${ARGN})
  list(LENGTH extra_args extra_count)
  # Empty list unless extra arguments are provided
  set(extra_depends )
  if (${extra_count} GREATER 0)
    # Set list to global variable named by the argument
    get_property(extra_depends GLOBAL PROPERTY ${ARGN})
  endif()
  blt_add_library(NAME Spheral_${package_name}
    HEADERS     ${${package_name}_headers}
    SOURCES     ${${package_name}_sources}
    DEPENDS_ON  ${spheral_blt_depends} ${${package_name}_ADDITIONAL_DEPENDS} ${SPHERAL_CXX_DEPENDS} ${extra_depends}
    OBJECT      TRUE)
  target_include_directories(Spheral_${package_name} SYSTEM PUBLIC ${SPHERAL_SUBMOD_INCLUDES})
  # Install the headers
  install(FILES ${${package_name}_headers}
    DESTINATION include/${package_name})

  if(ENABLE_CUDA)
    set_target_properties(Spheral_${package_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()

endfunction()

#----------------------------------------------------------------------------------------
#                                   spheral_add_cxx_library
#----------------------------------------------------------------------------------------
# -------------------------------------------
# VARIABLES THAT NEED TO BE PREVIOUSLY DEFINED
# -------------------------------------------
# spheral_blt_depends    : REQUIRED : List of external dependencies
# SPHERAL_CXX_DEPENDS    : REQUIRED : List of compiler dependencies
# <package_name>_headers : OPTIONAL : List of necessary headers to include
# <package_name>_sources : OPTIONAL : List of necessary source files to include
# SPHERAL_SUBMOD_DEPENDS : REQUIRED : List of submodule dependencies
# ----------------------
# INPUT-OUTPUT VARIABLES
# ----------------------
# <package_name>   : REQUIRED : Desired package name
# <_cxx_obj_list>  : REQUIRED : List of internal targets to include
# -----------------------
# OUTPUT VARIABLES TO USE - Made available implicitly after function call
# -----------------------
# <Spheral_package_name> : Exportable target for interal package name library
#----------------------------------------------------------------------------------------
function(spheral_add_cxx_library package_name _cxx_obj_list)
  # Assumes global variable spheral_blt_depends exists and is filled with external dependencies
  get_property(spheral_blt_depends GLOBAL PROPERTY spheral_blt_depends)
  # Assumes global variable spheral_cxx_depends exists and is filled with compiler dependencies
  get_property(SPHERAL_CXX_DEPENDS GLOBAL PROPERTY SPHERAL_CXX_DEPENDS)
  # For including files in submodules, currently unused
  get_property(SPHERAL_SUBMOD_INCLUDES GLOBAL PROPERTY SPHERAL_SUBMOD_INCLUDES)

  # set(EXTRA_CXX_DEPENDS )
  # # We can optionally specify the obj_libs_list and any additional dependencies
  # set(extra_args ${ARGN})
  # list(LENGTH extra_args extra_count)
  # if (${extra_count} GREATER 0)
  #   list(GET extra_args 0 spheral_blt_depends)
  # endif()
  # if (${extra_count} GREATER 1)
  #   list(GET extra_args 1 optional_arg)
  #   list(APPEND EXTRA_CXX_DEPENDS ${optional_arg})
  # endif()
  if(NOT ENABLE_SHARED)
    # Build static spheral C++ library
    blt_add_library(NAME Spheral_${package_name}
      HEADERS     ${${package_name}_headers}
      SOURCES     ${${package_name}_sources}
      DEPENDS_ON  ${_cxx_obj_list} ${spheral_blt_depends} ${SPHERAL_CXX_DEPENDS}
      SHARED      FALSE)
  else()
    # Build shared spheral C++ library
    blt_add_library(NAME Spheral_${package_name}
      HEADERS     ${${package_name}_headers}
      SOURCES     ${${package_name}_sources}
      DEPENDS_ON  ${_cxx_obj_list} ${spheral_blt_depends} ${SPHERAL_CXX_DEPENDS}
      SHARED      TRUE)
  endif()
  target_include_directories(Spheral_${package_name} SYSTEM PRIVATE ${SPHERAL_SUBMOD_INCLUDES})
  if(ENABLE_CUDA)
    set_target_properties(Spheral_${package_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()

  # Convert package name to lower-case for export target name
  string(TOLOWER ${package_name} lower_case_package)

  # Install Spheral C++ target and set it as an exportable CMake target
  install(TARGETS Spheral_${package_name}
    DESTINATION   lib
    EXPORT        spheral_${lower_case_package}-targets)
  # Export Spheral target
  install(EXPORT spheral_${lower_case_package}-targets DESTINATION lib/cmake)

  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  # TODO: Determine if this is still necessary
  set_target_properties(Spheral_${package_name} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  #${conduit_DIR}/lib;${axom_DIR}/lib;${boost_DIR}/lib;${hdf5_DIR}/lib;${zlib_DIR}/lib;${SPHERAL_ADDITIONAL_RPATHS}")
endfunction()


#-----------------------------------------------------------------------------------
# spheral_add_pybind11_library
#     - Generate the python friendly Spheral package lib
#
# Args:
#   package_name : *name* of spheral package to make into a library
#   INCLUDES     : optional, any additional include paths
#   SOURCES      : optional, any additional source files to compile into the library
#   DEPENDS      : optional, extra dependencies
# 
# Variables that must be set before calling spheral_add_obj_library:
#     spheral_depends
#         - List of targets the library depends on
#     spheral_blt_depends
#         - List of blt/libs the library depends on
#
#-----------------------------------------------------------------------------------

function(spheral_add_pybind11_library package_name)

  # Define our arguments
  set(options )
  set(oneValueArgs )
  set(multiValueArgs INCLUDES SOURCES DEPENDS)
  cmake_parse_arguments(${package_name} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  # message("** ${package_name}_INCLUDES: ${${package_name}_INCLUDES}")
  # message("** ${package_name}_SOURCES: ${${package_name}_SOURCES}")
  # message("** ${package_name}_DEPENDS: ${${package_name}_DEPENDS}")

  # List directories in which spheral .py files can be found.
  set(PYTHON_ENV 
      ${EXTRA_PYB11_SPHERAL_ENV_VARS}
      "${BUILDTIME_PYTHONENV_STR}:"
      "${SPHERAL_ROOT_DIR}/src/PYB11:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/${PYB11_MODULE_NAME}:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/polytope:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Distributed:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/OpenMP:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/CXXTypes:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Geometry:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/PolyClipper:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Silo:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/DataOutput:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/NodeList:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Field:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FieldList:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Kernel:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Neighbor:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Material:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FileIO:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/DataBase:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Boundary:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Physics:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Hydro:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ExternalForce:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Gravity:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Integrator:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Utilities:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/NodeGenerators:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/FieldOperations:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SPH:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/RK:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/CRKSPH:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ArtificialViscosity:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SVPH:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Mesh:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Damage:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/SolidMaterial:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/Strength:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/ArtificialConduction:"
      "${SPHERAL_ROOT_DIR}/src/PYB11/KernelIntegrator:"
      "${CMAKE_BINARY_DIR}/src/SimulationControl"
      )

  # Format list into a one line shell friendly format
  STRING(REPLACE ";" "<->" PYTHON_ENV_STR ${PYTHON_ENV})

  string(JOIN ":" PYTHON_ENV_STR ${PYTHON_ENV_STR} ${SPACK_PYTHONPATH})

  # Flag to tell PYB11Generator to use blt_add_library
  set(${package_name}_USE_BLT ON BOOL)
  # Get the TPL dependencies
  get_property(spheral_blt_depends GLOBAL PROPERTY spheral_blt_depends)
  get_property(spheral_tpl_includes GLOBAL PROPERTY spheral_tpl_includes)
  get_property(spheral_tpl_libraries GLOBAL PROPERTY spheral_tpl_libraries)

  set(MODULE_NAME Spheral${package_name})
  PYB11Generator_add_module(${package_name}
                            MODULE          ${MODULE_NAME}
                            SOURCE          ${package_name}_PYB11.py
                            DEPENDS         ${spheral_depends} ${spheral_blt_depends} ${${package_name}_DEPENDS} ${SPHERAL_CXX_DEPENDS} ${EXTRA_CXX_DEPENDS} Spheral_CXX
                            PYTHONPATH      ${PYTHON_ENV_STR}
                            INCLUDES        ${CMAKE_CURRENT_SOURCE_DIR} ${SPHERAL_INCLUDES} ${${package_name}_INCLUDES} ${spheral_tpl_includes} ${PYBIND11_ROOT_DIR}/include 
                            LINKS           ${spheral_tpl_libraries}
                            COMPILE_OPTIONS ${SPHERAL_PYB11_TARGET_FLAGS}
                            USE_BLT         ON
                            EXTRA_SOURCE    ${${package_name}_SOURCES}
                            INSTALL         OFF
                            )
  target_include_directories(${MODULE_NAME} SYSTEM PRIVATE ${SPHERAL_EXTERN_INCLUDES})
  target_compile_options(${MODULE_NAME} PRIVATE ${SPHERAL_PYB11_TARGET_FLAGS})

  install(TARGETS     ${MODULE_NAME}
          DESTINATION Spheral
          )

  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  set_target_properties(${MODULE_NAME} PROPERTIES
                        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${conduit_DIR}/lib;${axom_DIR}/lib;${boost_DIR}/lib;${python_DIR}/lib;${hdf5_DIR}/lib;${zlib_DIR}/lib;${SPHERAL_ADDITIONAL_RPATHS}"
                        )

endfunction()
