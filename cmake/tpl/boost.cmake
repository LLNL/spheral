if (NOT BOOST_HEADER_ONLY)
  set(${lib_name}_libs libboost_filesystem.so libboost_system.so)
endif()
