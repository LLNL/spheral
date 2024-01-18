if(APPLE)
  set(SHARED_EXT "dylib")
else()
  set(SHARED_EXT "so")
endif()

set(${lib_name}_libs libpolytope.${SHARED_EXT})

if(ENABLE_STATIC_TPL)
  string(REPLACE ".${SHARED_EXT}" ".a;" ${lib_name}_libs ${${lib_name}_libs})
endif()
