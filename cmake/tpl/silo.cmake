set(${lib_name}_libs libsiloh5.a)

if(APPLE)
  set(SHARED_EXT "dylib")
  set(${lib_name}_libs libsiloh5.${SHARED_EXT})
endif()
