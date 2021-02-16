set(POLYCLIPPER_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(POLYCLIPPER_DIST "PolyClipper-1.01.zip")
set(POLYCLIPPER_CACHE "${CACHE_DIR}/${POLYCLIPPER_DIST}")
set(POLYCLIPPER_URL "https://github.com/LLNL/PolyClipper/archive/v1.01.tar.gz")
set(POLYCLIPPER_MD5 "MD5=85b6b9e0dc47ed6fa4720f07cb1a7ff7")
set(POLYCLIPPER_DEST_DIR "${${lib_name}_DIR}/lib")

set(${lib_name}_libs libpolyclipper.a)

if(ENABLE_CXXONLY)
  set(POLYCLIPPER_ENABLE_CXXONLY On)
else()
  set(POLYCLIPPER_ENABLE_CXXONLY Off)
  list(APPEND POLYCLIPPER_DEPENDS python-install pip-modules ${spheral_py_depends})
endif()

if(${lib_name}_BUILD)

  if (EXISTS ${POLYCLIPPER_CACHE})
    set(POLYCLIPPER_URL ${POLYCLIPPER_CACHE})
  endif()

  ExternalProject_add(${lib_name}
    PREFIX ${POLYCLIPPER_PREFIX}
    URL ${POLYCLIPPER_URL}
    URL_HASH ${POLYCLIPPER_MD5}
    DOWNLOAD_DIR ${CACHE_DIR}
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_INSTALL_PREFIX=${${lib_name}_DIR}
               -DPOLYCLIPPER_BLT_DIR=${CMAKE_SOURCE_DIR}/cmake/blt
               -DENABLE_CXXONLY=${POLYCLIPPER_ENABLE_CXXONLY}
               -DPYTHON_EXE=${PYTHON_EXE}
               -DLOOKUP_PYBIND11_INCLUDE_PATH=On
               -DPOLYCLIPPER_PYTHON_INSTALL=${PYTHON_SITE_PACKAGE_DIR}
               -DENABLE_DOCS=Off
               DEPENDS ${POLYCLIPPER_DEPENDS}
    LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
    LOG_CONFIGURE ${OUT_PROTOCOL_EP}
    LOG_BUILD ${OUT_PROTOCOL_EP}
    LOG_INSTALL ${OUT_PROTOCOL_EP}
  )

  # if(NOT ENABLE_CXXONLY AND NOT ENABLE_STATIC_CXXONLY)
  #   install(
  #     FILES ${${lib_name}_DIR}/lib/python2.7/site-packages/polyclipper/polyclipper.so
  #     DESTINATION ${PYTHON_SITE_PACKAGE_DIR}
  #     )
  # endif()
endif()

