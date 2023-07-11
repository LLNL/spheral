# Ignore -Wterminate warning that is thrown from DBC.hh VERIFY2

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set_source_files_properties(SpheralSilo.cc PROPERTIES COMPILE_FLAGS -Wno-terminate)
endif()

spheral_add_pybind11_library(Silo
                             INCLUDES ${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Silo)
