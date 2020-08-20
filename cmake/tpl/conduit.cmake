set(CONDUIT_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(CONDUIT_DIST "conduit-v0.5.1-src-with-blt.tar.gz")
set(CONDUIT_URL "https://github.com/LLNL/conduit/releases/download/v0.5.1/${CONDUIT_DIST}")
set(CONDUIT_CACHE "${CACHE_DIR}/${CONDUIT_DIST}")

set(${lib_name}_libs )

if(${lib_name}_BUILD)

  if (EXISTS ${CONDUIT_CACHE})
    set(CONDUIT_URL ${CONDUIT_CACHE})
  endif()

  ExternalProject_add(${lib_name}
    PREFIX ${CONDUIT_PREFIX}
    URL ${CONDUIT_URL}
    DOWNLOAD_DIR ${CACHE_DIR}
    CMAKE_ARGS ../conduit/src/ -DCMAKE_BUILD_TYPE=Release
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_C_FLAGS=-fPIC
               -DENABLE_TESTS=Off
               -DHDF5_DIR=${HDF5_INSTALL_DIR}
               -DCMAKE_INSTALL_PREFIX=${${lib_name}_DIR}

    DEPENDS hdf5

    LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
    LOG_CONFIGURE ${OUT_PROTOCOL_EP}
    LOG_BUILD ${OUT_PROTOCOL_EP}
    LOG_INSTALL ${OUT_PROTOCOL_EP}
  )

endif()
set(CONDUIT_INSTALL_DIR ${${lib_name}_DIR} PARENT_SCOPE)
