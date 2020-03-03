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


