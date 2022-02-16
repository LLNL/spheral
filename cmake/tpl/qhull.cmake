set(QHULL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(QHULL_URL "https://github.com/qhull/qhull/archive/2019.1.tar.gz")
set(QHULL_CACHE "${CACHE_DIR}/2019.1.tar.gz")
set(QHULL_SRC_DIR ${QHULL_PREFIX}/src/qhull/src)

# Setting this to just the release library until we support TPL debug builds on LC
set(${lib_name}_libs libqhull.so)
if(ENABLE_STATIC_TPL)
  set(${lib_name}_libs libqhullstatic.a)
endif()
set(QHULL_BUILD_TYPE Release)
# if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#   set(${lib_name}_libs libqhullstatic_d.a)
# else()
#   set(${lib_name}_libs libqhullstatic.a)
# endif()
# set(QHULL_BUILD_TYPE ${CMAKE_BUILD_TYPE})

if(${lib_name}_BUILD)

  if (EXISTS ${QHULL_CACHE})
    set(QHULL_URL ${QHULL_CACHE})
  endif()

  ExternalProject_add(${lib_name}
    PREFIX ${QHULL_PREFIX}
    PATCH_COMMAND patch -t ${QHULL_SRC_DIR}/libqhull/qhull_a.h ${PATCH_DIR}/qhull-2015.2-qhull_a.h-patch &&
                  patch -t ${QHULL_SRC_DIR}/libqhull_r/qhull_ra.h ${PATCH_DIR}/qhull-2015.2-qhull_ra.h-patch
    URL ${QHULL_URL}
    URL_HASH "MD5=${QHULL_MD5}"
    DOWNLOAD_DIR ${CACHE_DIR}
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_C_FLAGS=-fPIC
               -DBUILD_SHARED_LIBS=OFF
               -DCMAKE_INSTALL_PREFIX=${${lib_name}_DIR}
               -DCMAKE_BUILD_TYPE=${QHULL_BUILD_TYPE}

    LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
    LOG_CONFIGURE ${OUT_PROTOCOL_EP}
    LOG_BUILD ${OUT_PROTOCOL_EP}
    LOG_INSTALL ${OUT_PROTOCOL_EP}
  )

endif()
