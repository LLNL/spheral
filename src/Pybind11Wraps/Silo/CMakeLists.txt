# Ignore -Wterminate warning that is thrown from DBC.hh VERIFY2

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
if (NOT ENABLE_CUDA)
  set_source_files_properties(SpheralSilo.cc PROPERTIES COMPILE_FLAGS -Wno-terminate)
endif()
endif()

set(Silo_ADDITIONAL_INCLUDES ${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Silo)
spheral_add_pybind11_library(Silo)
