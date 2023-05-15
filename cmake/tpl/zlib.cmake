set(${lib_name}_libs libz.so)

if(ENABLE_STATIC_TPL)
  string(REPLACE ".so" ".a;" ${lib_name}_libs ${${lib_name}_libs})
endif()
