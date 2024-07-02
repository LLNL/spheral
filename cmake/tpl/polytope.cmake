set(${lib_name}_libs libpolytope.a)

if(APPLE)
  set(SHARED_EXT "dylib")
  set(${lib_name}_libs libpolytope.${SHARED_EXT})
endif()
