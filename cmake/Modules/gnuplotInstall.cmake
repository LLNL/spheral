set(GNUPLOT_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(GNUPLOT_DIST gnuplot-py-1.8.tar.gz)
set(GNUPLOT_URL http://downloads.sourceforge.net/gnuplot-py/gnuplot-py-1.8.tar.gz)
set(GNUPLOT_HASH )
set(GNUPLOT_CACHE ${CACHE_DIR}/${GNUPLOT_DIST})

if(EXISTS ${GNUPLOT_CACHE})
  set(GNUPLOT_URL ${GNUPLOT_CACHE})
endif()

ExternalProject_add(${lib_name}
  URL ${GNUPLOT_URL} 
  PREFIX ${GNUPLOT_PREFIX}
  DOWNLOAD_DIR ${CACHE_DIR}
  CONFIGURE_COMMAND ${PYTHON_EXE} ${GNUPLOT_PREFIX}/src/gnuplot/setup.py install
  BUILD_COMMAND sleep 1
  INSTALL_COMMAND sleep 1
  DEPENDS python-install numpy-stl
  LOG_DOWNLOAD ${OUT_PROTOCOL_EP}
  LOG_CONFIGURE ${OUT_PROTOCOL_EP}
)
