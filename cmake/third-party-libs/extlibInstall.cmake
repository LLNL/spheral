set(${lib_name}_libs libmylib.a)

if(${lib_name}_BUILD)
ExternalProject_Add(${lib_name}
  SOURCE_DIR "/usr/workspace/wsrzd/davis291/CppExamples/CMAKE_Examples/mylib"
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${${lib_name}_DIR}
  )
endif()
