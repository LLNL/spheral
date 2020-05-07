function(spheral_add_cxx_library package_name)
  if(ENABLE_STATIC_CXXONLY)
    blt_add_library(NAME        Spheral_${package_name}
                    HEADERS     ${${package_name}_headers}
                    SOURCES     ${${package_name}_sources}
                    DEPENDS_ON  -Wl,--start-group ${spheral_blt_depends} -Wl,--end-group
                    SHARED      FALSE
                    )
  else()
    blt_add_library(NAME        Spheral_${package_name}
                    HEADERS     ${${package_name}_headers}
                    SOURCES     ${${package_name}_sources}
                    DEPENDS_ON  -Wl,--start-group  ${spheral_blt_depends} -Wl,--end-group
                    SHARED      TRUE
                    )
  endif()

  if(spheral_depends)
    add_dependencies(Spheral_${package_name} ${spheral_depends})
  endif()

  install(TARGETS             Spheral_${package_name}
          EXPORT              ${PROJECT_NAME}-targets
          DESTINATION         lib
          INCLUDES DESTINATION include
          )

  install(FILES       ${${package_name}_headers}
          DESTINATION include/${package_name}
          )

#  set_target_properties(Spheral_${package_name} PROPERTIES
#    INSTALL_RPATH_USE_LINK_PATH TRUE)

endfunction()



function(spheral_add_pybind11_library package_name)
  include(${CMAKE_MODULE_PATH}/spheral/PYB11Generator.cmake)

  set(PYB11_MODULE_NAME ${package_name})
  PYB11_GENERATE_BINDINGS()

  set(MODULE_NAME Spheral${PYB11_MODULE_NAME})
  set(GENERATED_SOURCE Spheral${PYB11_GENERATED_SOURCE})

  blt_add_library(
    NAME         ${MODULE_NAME}
    SOURCES      ${GENERATED_SOURCE} ${${package_name}_ADDITIONAL_SOURCES}
    DEPENDS_ON   -Wl,--start-group ${SPHERAL_PYTHON_DEPENDS} -Wl,--end-group ${${package_name}_ADDITIONAL_DEPENDS} ${spheral_blt_depends}
    INCLUDES     ${${package_name}_ADDITIONAL_INCLUDES}
    OUTPUT_NAME  ${MODULE_NAME}
    CLEAR_PREFIX TRUE
    SHARED       TRUE
    )
  add_dependencies(${MODULE_NAME} ${spheral_py_depends} ${spheral_depends})

  install(
    FILES ${CMAKE_BINARY_DIR}/lib/$<TARGET_FILE_NAME:${MODULE_NAME}>
    DESTINATION .
    )

  #  set_target_properties(${MODULE_NAME} PROPERTIES
  #    INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}
  #    INSTALL_RPATH_USE_LINK_PATH TRUE)

endfunction()
