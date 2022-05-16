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


