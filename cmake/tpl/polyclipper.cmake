set(POLYCLIPPER_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(POLYCLIPPER_DIST "PolyClipper-v.1.2.3.zip")
set(POLYCLIPPER_CACHE "${CACHE_DIR}/${POLYCLIPPER_DIST}")
set(POLYCLIPPER_URL "https://github.com/LLNL/PolyClipper/archive/refs/tags/v1.2.3.zip")
set(POLYCLIPPER_DEST_DIR "${${lib_name}_DIR}/lib")

set(${lib_name}_libs )

# We currently build our own Python bindings in Spheral since we use our own Vector types,
# so for simplicity I'm suspending building the built-in PolyClipper Python bindings.
set(POLYCLIPPER_ENABLE_CXXONLY On)
# if(ENABLE_CXXONLY)
#   set(POLYCLIPPER_ENABLE_CXXONLY On)
# else()
#   set(POLYCLIPPER_ENABLE_CXXONLY Off)
#   list(APPEND POLYCLIPPER_DEPENDS python-install ${spheral_py_depends})
# endif()

if(${lib_name}_BUILD)

  if (EXISTS ${POLYCLIPPER_CACHE})
    set(POLYCLIPPER_URL ${POLYCLIPPER_CACHE})
  endif()

  set(PYBIND11_INCLUDE_PATH "${PYTHON_INSTALL_DIR}/include/python2.7 ${PYBIND11_INSTALL_DIR}/include")

  ExternalProject_add(${lib_name}
    PREFIX ${POLYCLIPPER_PREFIX}
    URL ${POLYCLIPPER_URL}
    URL_HASH "MD5=${POLYCLIPPER_MD5}"
    DOWNLOAD_DIR ${CACHE_DIR}
    DEPENDS ${pybind11_build_dep} ${pip-modules_build_dep}
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_INSTALL_PREFIX=${${lib_name}_DIR}
               -DPOLYCLIPPER_BLT_DIR=${SPHERAL_BLT_DIR}
               -DENABLE_CXXONLY=${POLYCLIPPER_ENABLE_CXXONLY}
               -DENABLE_OPENMP=${ENABLE_OPENMP}
               -DENABLE_MPI=Off
               -DPYTHON_EXE=${PYTHON_EXE}
               #-DLOOKUP_PYBIND11_INCLUDE_PATH=On
               -DPYBIND11_INCLUDE_PATH=${PYBIND11_INCLUDE_PATH}
               -DPYB11GEN_PATH=${PYTHON_SITE_PACKAGE_DIR}
               -DPOLYCLIPPER_PYTHON_INSTALL=${${lib_name}_DIR}
               -DENABLE_DOCS=Off
               DEPENDS ${POLYCLIPPER_DEPENDS}
    LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
    LOG_CONFIGURE ${OUT_PROTOCOL_EP}
    LOG_BUILD ${OUT_PROTOCOL_EP}
    LOG_INSTALL ${OUT_PROTOCOL_EP}
  )
endif()

