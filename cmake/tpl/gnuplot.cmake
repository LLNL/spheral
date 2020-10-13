set(GNUPLOT_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(GNUPLOT_DIST gnuplot-py-1.8.tar.gz)
set(GNUPLOT_URL http://downloads.sourceforge.net/gnuplot-py/gnuplot-py-1.8.tar.gz)
set(GNUPLOT_HASH )
set(GNUPLOT_CACHE ${CACHE_DIR}/${GNUPLOT_DIST})

set(GNUPLOT_DOWNLOAD_CMD wget ${GNUPLOT_URL} --no-check-certificate -O ${CACHE_DIR}/gnuplot-py-1.8.tar.gz)

if(python_BUILD)
  if(EXISTS ${GNUPLOT_CACHE})
    set(GNUPLOT_DOWNLOAD_CMD sleep 1)
  endif()

  add_custom_target(${lib_name}-unpack
    COMMAND ${GNUPLOT_DOWNLOAD_CMD} > gnuplot-download.log
    COMMAND tar -xvf ${CACHE_DIR}/gnuplot-py-1.8.tar.gz > gnuplot-unpack.log
    COMMENT "Upacking gnuplot"
  )

  add_custom_target(${lib_name}
    COMMAND ${PYTHON_EXE} setup.py install > gnuplot-install.log
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gnuplot-py-1.8
    COMMENT "Installing gnuplot"
    DEPENDS python-install pip-modules ${lib_name}-unpack
  )
endif()
