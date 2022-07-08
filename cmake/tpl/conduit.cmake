list(APPEND ${lib_name}_INCLUDES $<BUILD_INTERFACE:${${lib_name}_DIR}/include/${lib_name}>)

set(${lib_name}_libs 
    libconduit.a
    libconduit_blueprint.a
    libconduit_relay.a
   )

if(ENABLE_MPI)
  list(APPEND ${lib_name}_libs
    libconduit_blueprint_mpi.a
    libconduit_relay_mpi.a
    libconduit_relay_mpi_io.a)
endif()

#if(ENABLE_STATIC_TPL)
#  string(REPLACE ".so" ".a;" ${lib_name}_libs ${${lib_name}_libs})
#endif()
