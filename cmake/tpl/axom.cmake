set(AXOM_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(AXOM_DIST "Axom-v0.3.3.tar.gz")
set(AXOM_URL "https://github.com/LLNL/axom/releases/download/v0.3.3/${AXOM_DIST}")
set(AXOM_CACHE "${CACHE_DIR}/${AXOM_DIST}")

set(${lib_name}_libs 
    libaxom.a 
   )
