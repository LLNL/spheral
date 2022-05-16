set(OPENSUBDIV_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(OPENSUBDIV_DIST "v3_3_0.zip")
set(OPENSUBDIV_URL "https://github.com/PixarAnimationStudios/OpenSubdiv/archive/${OPENSUBDIV_DIST}")
set(OPENSUBDIV_CACHE "${CACHE_DIR}/${OPENSUBDIV_DIST}")
set(OPENSUBDIV_SRC_DIR ${OPENSUBDIV_PREFIX}/src/opensubdiv/src)

set(${lib_name}_libs libosdCPU.a)

