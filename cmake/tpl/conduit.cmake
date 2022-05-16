set(CONDUIT_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(CONDUIT_DIST "conduit-v0.5.1-src-with-blt.tar.gz")
set(CONDUIT_URL "https://github.com/LLNL/conduit/releases/download/v0.5.1/${CONDUIT_DIST}")
set(CONDUIT_CACHE "${CACHE_DIR}/${CONDUIT_DIST}")

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

set(CONDUIT_INSTALL_DIR ${${lib_name}_DIR} PARENT_SCOPE)
