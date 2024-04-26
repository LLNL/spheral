set(${lib_name}_libs libsiloh5.a)

if(APPLE)
  set(SHARED_EXT "dylib")
  set(${lib_name}_libs libsiloh5.${SHARED_EXT})
endif()

if(ENABLE_STATIC_TPL)
  string(REPLACE ".${SHARED_EXT}" ".a;" ${lib_name}_libs ${${lib_name}_libs})
endif()
