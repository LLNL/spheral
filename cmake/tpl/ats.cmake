set(ATS_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${lib_name})
set(ATS_DIST ats-5.2.tar.gz)
set(ATS_URL "https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/ats/${ATS_DIST}")
set(ATS_HASH )
set(ATS_CACHE ${CACHE_DIR}/${ATS_DIST})

set(ATS_DOWNLOAD_CMD wget ${ATS_URL} --no-check-certificate -O ${CACHE_DIR}/${ATS_DIST})

if(python_BUILD)
  if(EXISTS ${ATS_CACHE})
    set(ATS_DOWNLOAD_CMD sleep 1)
  endif()

  add_custom_target(${lib_name}-unpack
    COMMAND ${ATS_DOWNLOAD_CMD} > ats-download.log
    COMMAND tar -zxvf ${CACHE_DIR}/${ATS_DIST} > ats-unpack.log
    COMMENT "Upacking ats"
  )

  add_custom_target(${lib_name}
    COMMAND ${PYTHON_EXE} setup.py install > ats-install.log
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ats-5.2
    COMMENT "Installing ats"
    DEPENDS python-install pip-modules ${lib_name}-unpack
  )

  list(APPEND spheral_py_depends ${lib_name})
endif()
