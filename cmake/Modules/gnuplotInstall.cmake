set(GNUPLOT_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(GNUPLOT_DIST gnuplot-py-1.8.tar.gz)
set(GNUPLOT_URL http://downloads.sourceforge.net/gnuplot-py/gnuplot-py-1.8.tar.gz)
set(GNUPLOT_HASH )
set(GNUPLOT_CACHE ${CACHE_DIR}/${GNUPLOT_DIST})

set(GNUPLOT_DOWNLOAD_CMD wget ${GNUPLOT_URL} -O ${CACHE_DIR}/gnuplot-py-1.8.tar.gz)

if(EXISTS ${GNUPLOT_CACHE})
  set(GNUPLOT_DOWNLOAD_CMD sleep 1)
endif()

add_custom_target(${lib_name}-unpack
  COMMAND ${GNUPLOT_DOWNLOAD_CMD}
  COMMAND tar -xvf ${CACHE_DIR}/gnuplot-py-1.8.tar.gz
  COMMENT "Upacking gnuplot"
)

add_custom_target(${lib_name}
  COMMAND ${PYTHON_EXE} setup.py install
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gnuplot-py-1.8
  COMMENT "Installing gnuplot"
  DEPENDS python-install pip-modules ${lib_name}-unpack
)
