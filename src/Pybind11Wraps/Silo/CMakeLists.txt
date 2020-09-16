# Ignore -Wterminate warning that is thrown from DBC.hh VERIFY2
set_source_files_properties(SpheralSilo.cc PROPERTIES COMPILE_FLAGS -Wno-terminate)

set(Silo_ADDITIONAL_INCLUDES ${SPHERAL_ROOT_DIR}/src/Pybind11Wraps/Silo)
spheral_add_pybind11_library(Silo)
