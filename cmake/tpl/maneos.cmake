set(MANEOS_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(MANEOS_DIST "M-ANEOS-v1.0.tar.gz")
set(MANEOS_CACHE "${CACHE_DIR}/${MANEOS_DIST}")
set(MANEOS_URL "https://github.com/isale-code/M-ANEOS/releases/download/v1.0beta/${MANEOS_DIST}")
set(MANEOS_SRC_DIR "${MANEOS_PREFIX}/src/maneos/src")
set(MANEOS_DEST_DIR "${SPHERAL_INSTALL_DIR}/maneos/lib")

#set(${lib_name}_INCLUDES maneos.h)
set(${lib_name}_libs libaneos.a)

if(${lib_name}_BUILD)

  if (EXISTS ${MANEOS_CACHE})
    set(MANEOS_URL ${MANEOS_CACHE})
  endif()

  ExternalProject_add(${lib_name}
    PREFIX ${MANEOS_PREFIX}
    URL ${MANEOS_URL} 
    DOWNLOAD_DIR ${CACHE_DIR}
    PATCH_COMMAND patch -t ${MANEOS_SRC_DIR}/Makefile ${PATCH_DIR}/MANEOS-v1.0_Makefile.patch
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make FORT=${CMAKE_Fortran_COMPILER} -C ${MANEOS_SRC_DIR}
    BUILD_IN_SOURCE 1
    INSTALL_DIR ${MANEOS_DEST_DIR}
    INSTALL_COMMAND cp -f ${MANEOS_SRC_DIR}/libaneos.a ${MANEOS_DEST_DIR}
    LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
    LOG_CONFIGURE ${OUT_PROTOCOL_EP}
    LOG_BUILD ${OUT_PROTOCOL_EP}
    LOG_INSTALL ${OUT_PROTOCOL_EP}
  )

endif()

install(
  FILES ${MANEOS_SRC_DIR}/libaneos.a
  TYPE LIB)
