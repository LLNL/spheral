set(POLYCLIPPER_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(POLYCLIPPER_DIST "PolyClipper-1.01.zip")
set(POLYCLIPPER_CACHE "${CACHE_DIR}/${POLYCLIPPER_DIST}")
set(POLYCLIPPER_URL "https://github.com/LLNL/PolyClipper/archive/v1.01.tar.gz")
set(POLYCLIPPER_MD5 "MD5=2a249c4d0350f4d280c9effacafd4a9a")
set(POLYCLIPPER_DEST_DIR "${${lib_name}_DIR}/lib")





set(POLYCLIPPER_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(POLYCLIPPER_CACHE "${CACHE_DIR}/0.6.2.tar.gz")
set(POLYCLIPPER_URL "https://github.com/pbtoast/polyclipper/archive/0.6.2.tar.gz")
set(POLYCLIPPER_MD5 "MD5=c8cad6ee08fbafaf074359d28486cf2c")
set(POLTOPE_SRC_DIR ${POLYCLIPPER_PREFIX}/src/boost)

set(${lib_name}_libs libpolyclipper.a)

set(POLYCLIPPER_DEPENDS boost)
set(POLYCLIPPER_USE_PYTHON On)

if(ENABLE_CXXONLY)
  set(POLYCLIPPER_USE_PYTHON Off)
else()
  list(APPEND POLYCLIPPER_DEPENDS python-install pip-modules ${spheral_py_depends})
endif()

if(${lib_name}_BUILD)

  if (EXISTS ${POLYCLIPPER_CACHE})
    set(POLYCLIPPER_URL ${POLYCLIPPER_CACHE})
  endif()

  ExternalProject_add(${lib_name}
    PREFIX ${POLYCLIPPER_PREFIX}
    PATCH_COMMAND patch -t ${POLYCLIPPER_PREFIX}/src/polyclipper/src/PYB11/CMakeLists.txt  ${PATCH_DIR}/polyclipper-PYB11-CMakeLists.patch

    URL ${POLYCLIPPER_URL}
    URL_HASH ${POLYCLIPPER_MD5}
    DOWNLOAD_DIR ${CACHE_DIR}
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_INSTALL_PREFIX=${${lib_name}_DIR} 
               -DPYBIND11_INCLUDE_DIRS=${PYBIND11_INSTALL_DIR}/include
               -DUSE_PYTHON=${POLYCLIPPER_USE_PYTHON}
               -DPYTHON_EXE=${PYTHON_EXE}
               -DBoost_INCLUDE_DIR=${BOOST_INSTALL_DIR}/include
               -DTESTING=Off
               DEPENDS ${POLYCLIPPER_DEPENDS}

    LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
    LOG_CONFIGURE ${OUT_PROTOCOL_EP}
    LOG_BUILD ${OUT_PROTOCOL_EP}
    LOG_INSTALL ${OUT_PROTOCOL_EP}
  )

  if(NOT ENABLE_CXXONLY AND NOT ENABLE_STATIC_CXXONLY)
    install(
      FILES ${${lib_name}_DIR}/lib/python2.7/site-packages/polyclipper/polyclipper.so
      DESTINATION ${PYTHON_SITE_PACKAGE_DIR}
      )
  endif()
endif()

