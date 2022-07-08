list(APPEND ${lib_name}_INCLUDES $<BUILD_INTERFACE:${${lib_name}_DIR}/include/${lib_name}>)

set(${lib_name}_libs 
    libconduit.so
    libconduit_blueprint.so
    libconduit_relay.so
   )

if(ENABLE_MPI)
  list(APPEND ${lib_name}_libs
    libconduit_blueprint_mpi.so
    libconduit_relay_mpi.so
    libconduit_relay_mpi_io.so)
endif()

if(ENABLE_STATIC_TPL)
  string(REPLACE ".so" ".a;" ${lib_name}_libs ${${lib_name}_libs})
endif()
