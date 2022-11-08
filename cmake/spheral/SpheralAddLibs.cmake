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
                  DEPENDS_ON  ${spheral_blt_depends} ${spheral_blt_cxx_depends} blt_python ${${package_name}_ADDITIONAL_DEPENDS} ${SPHERAL_CXX_DEPENDS}
                  OBJECT TRUE
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
                        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${conduit_DIR}/lib;${axom_DIR}/lib;${boost_DIR}/lib"
                        )
endfunction()


#-----------------------------------------------------------------------------------
# spheral_add_pybind11_library
#     - Generate the python friendly Spheral package lib
#
# package_name : *name* of spheral package to make into a library
#
# Variables that must be set before calling spheral_add_obj_library:
#     <package_Name>_ADDITIONAL_INCLUDES
#         - List of addition includes needed by a given package
#     <package_name>_ADDITIONAL_SOURCE
#         - List of additional sources to build library with
#     spheral_depends
#         - List of targets the library depends on
#     spheral_blt_depends
#         - List of blt/libs the library depends on
#
#-----------------------------------------------------------------------------------

function(spheral_add_pybind11_library package_name)
  include(${CMAKE_MODULE_PATH}/spheral/PYB11Generator.cmake)

  # Generate the pybind11 C++ source file
  PYB11_GENERATE_BINDINGS(${package_name})

  # Build python friendly spheral lib
  set(MODULE_NAME Spheral${package_name})
  blt_add_library(NAME         ${MODULE_NAME}
                  SOURCES      ${PYB11_GENERATED_SOURCE} ${${package_name}_ADDITIONAL_SOURCES}
                  DEPENDS_ON   Spheral_CXX ${spheral_blt_depends} ${spheral_blt_py_depends} ${${package_name}_ADDITIONAL_DEPENDS}
                  INCLUDES     ${${package_name}_ADDITIONAL_INCLUDES}
                  OUTPUT_NAME  ${MODULE_NAME}
                  CLEAR_PREFIX TRUE
                  SHARED       TRUE
                  )

  target_compile_options(${MODULE_NAME} PRIVATE ${SPHERAL_PYB11_TARGET_FLAGS})

  install(TARGETS     ${MODULE_NAME}
          DESTINATION Spheral
          )

  # Set the r-path of the C++ lib such that it is independent of the build dir when installed
  set_target_properties(${MODULE_NAME} PROPERTIES
                        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${conduit_DIR}/lib;${axom_DIR}/lib;${boost_DIR}/lib;${python_DIR}/lib"
                        )

endfunction()
