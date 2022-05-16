set(SILO_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(SILO_CACHE "${CACHE_DIR}/silo-4.10.2-bsd.tgz")
set(SILO_URL "https://wci.llnl.gov/sites/wci/files/2021-01/silo-4.10.2-bsd.tgz")
set(SILO_SRC_DIR ${SILO_PREFIX}/src/silo)

#set(${lib_name}_INCLUDES silo.h)
set(${lib_name}_libs libsiloh5.a)


