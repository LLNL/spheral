#-----------------------------------------------------------------------------------
# spheral_add_obj_library
#     - Generates the blt libraries for a given spheral C++ package
#
# package_name : *name* of spheral package to make into a library
# obj_libs_list : *name* of a CMake list which we will add the library object files to
#
# Variables that must be set before calling spheral_add_obj_library:
#     ENABLE_STATIC_CXXONLY : Default False
#         - If set to true all libs will be made static
#     <package_Name>_headers
#         - List of headers to be installed
#     <package_name>_sources
#         - List of sources to build library with
#     spheral_depends
#         - List of targets the library depends on
#     spheral_blt_depends
#         - List of blt/libs the library depends on
#
#-----------------------------------------------------------------------------------
function(spheral_add_obj_library
         package_name)

  set(obj_libs_list SPHERAL_OBJ_LIBS)

  # We can optionally specify the obj_libs_list
  set(extra_args ${ARGN})
  list(LENGTH extra_args extra_count)
  if (${extra_count} GREATER 0)
    list(GET extra_args 0 obj_libs_list)
  endif()

  blt_add_library(NAME        Spheral_${package_name}
                  HEADERS     ${${package_name}_headers}
                  SOURCES     ${${package_name}_sources}
                  DEPENDS_ON  ${spheral_blt_depends} ${spheral_blt_cxx_depends} ${${package_name}_ADDITIONAL_DEPENDS} ${SPHERAL_CXX_DEPENDS}
                  OBJECT TRUE
                  )
  target_include_directories(Spheral_${package_name} PRIVATE ${SPHERAL_INCLUDES})
  target_include_directories(Spheral_${package_name} SYSTEM PRIVATE ${SPHERAL_EXTERN_INCLUDES})

  # Install the headers
  install(FILES       ${${package_name}_headers}
          DESTINATION include/${package_name}
          )

  if(ENABLE_CUDA)
    set_target_properties(Spheral_${package_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()

  # Add this to the obj_libs_list list
  get_property(${obj_libs_list} GLOBAL PROPERTY ${obj_libs_list})
  list(APPEND ${obj_libs_list} Spheral_${package_name})
  set_property(GLOBAL PROPERTY ${obj_libs_list} "${${obj_libs_list}}")

endfunction()

#-----------------------------------------------------------------------------------
# spheral_add_cxx_library
#     - same interface as spheral_add_obj_library
#-----------------------------------------------------------------------------------
function(spheral_add_cxx_library
         package_name)

  set(obj_libs_list SPHERAL_OBJ_LIBS)
  set(EXTRA_CXX_DEPENDS )

  # We can optionally specify the obj_libs_list and any additional dependencies
  set(extra_args ${ARGN})
  list(LENGTH extra_args extra_count)
  if (${extra_count} GREATER 0)
    list(GET extra_args 0 obj_libs_list)
  endif()
  if (${extra_count} GREATER 1)
    list(GET extra_args 1 optional_arg)
    list(APPEND EXTRA_CXX_DEPENDS ${optional_arg})
  endif()
  get_property(${obj_libs_list} GLOBAL PROPERTY ${obj_libs_list})

  if(NOT ENABLE_SHARED)
    # Build static spheral C++ library
    blt_add_library(NAME        Spheral_${package_name}
                    HEADERS     ${${package_name}_headers}
                    SOURCES     ${${package_name}_sources}
                    DEPENDS_ON  ${${obj_libs_list}} ${SPHERAL_CXX_DEPENDS} ${EXTRA_CXX_DEPENDS}
                    SHARED      FALSE
                    )
  else()
    # Build shared spheral C++ library
    blt_add_library(NAME        Spheral_${package_name}
                    HEADERS     ${${package_name}_headers}
                    SOURCES     ${${package_name}_sources}
                    DEPENDS_ON  ${${obj_libs_list}} ${SPHERAL_CXX_DEPENDS} ${EXTRA_CXX_DEPENDS}
                    SHARED      TRUE
                    )
  endif()

  get_target_property(_INTERFACE_INCLUDE_DIRECTORIES Spheral_${package_name} INTERFACE_INCLUDE_DIRECTORIES)
  set_target_properties(Spheral_${package_name} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${_INTERFACE_INCLUDE_DIRECTORIES}>")

  get_target_property(_LINK_LIBRARIES Spheral_${package_name} LINK_LIBRARIES)
  LIST(REMOVE_DUPLICATES _LINK_LIBRARIES)
  set_target_properties(Spheral_${package_name} PROPERTIES LINK_LIBRARIES "${_LINK_LIBRARIES};${_LINK_LIBRARIES}") 

  if(ENABLE_CUDA)
    set_target_properties(Spheral_${package_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()

  # Install Spheral C++ target and set it as an exportable CMake target
  install(TARGETS             Spheral_${package_name}
          DESTINATION         lib
          EXPORT              ${PROJECT_NAME}-targets
          )

  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  set_target_properties(Spheral_${package_name} PROPERTIES
                        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${conduit_DIR}/lib;${axom_DIR}/lib;${boost_DIR}/lib;${hdf5_DIR}/lib"
                        )
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

  # Get the TPL dependencies
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
                            )
  target_include_directories(${package_name} SYSTEM PRIVATE ${SPHERAL_EXTERN_INCLUDES})
  target_compile_options(${package_name} PRIVATE ${SPHERAL_PYB11_TARGET_FLAGS})

  install(TARGETS     ${package_name}
          DESTINATION Spheral
          )

  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  set_target_properties(${package_name} PROPERTIES
                        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${conduit_DIR}/lib;${axom_DIR}/lib;${boost_DIR}/lib;${python_DIR}/lib;${hdf5_DIR}/lib"
                        )

endfunction()
