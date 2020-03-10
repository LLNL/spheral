function(spheral_add_cxx_library package_name)
  if(ENABLE_STATIC_CXXONLY)
    blt_add_library(NAME        Spheral_${package_name}
                    HEADERS     ${${package_name}_headers}
                    SOURCES     ${${package_name}_sources}
                    DEPENDS_ON  ${spheral_depends}
                    SHARED      FALSE
                    )
  else()
    blt_add_library(NAME        Spheral_${package_name}
                    HEADERS     ${${package_name}_headers}
                    SOURCES     ${${package_name}_sources}
                    DEPENDS_ON  ${spheral_depends}
                    SHARED      TRUE
                    )
  endif()
endfunction()

function(spheral_add_pybind11_library package_name)
  set(PYB11_MODULE_NAME ${package_name})
  include(UsePYB11Generator)
    PYB11_GENERATE_BINDINGS()

  set(MODULE_NAME Spheral${PYB11_MODULE_NAME})
  set(GENERATED_SOURCE Spheral${PYB11_GENERATED_SOURCE})

  blt_add_library(
    NAME         ${MODULE_NAME}
    SOURCES      ${GENERATED_SOURCE} ${${package_name}_ADDITIONAL_SOURCES}
    DEPENDS_ON   ${SPHERAL_PYTHON_DEPENDS} ${${package_name}_ADDITIONAL_DEPENDS} ${spheral_depends}
    INCLUDES     ${${package_name}_ADDITIONAL_INCLUDES}
    OUTPUT_NAME  ${MODULE_NAME}
    CLEAR_PREFIX TRUE
    SHARED       TRUE
    )

  install(
    FILES ${CMAKE_BINARY_DIR}/lib/$<TARGET_FILE_NAME:${MODULE_NAME}>
    DESTINATION .
    )
endfunction()
