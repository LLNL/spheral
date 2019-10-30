set(FileIO_inst
  )

set(FileIO_sources
  #DbFileIO.cc
  FileIO.cc
  #FlatFileIO.cc
  #HDF5IO.cc
  #HDF5Traits.cc
  #HDF5Types.cc
  #PyFileIO.cc
  #SiloFileIO.cc
  #vectorstringUtilities.cc
  )

instantiate(FileIO_inst FileIO_sources)

set(FileIO_headers
  DbFileIO.hh
  FileIO.hh
  FileIOInline.hh
  FlatFileIO.hh
  HDF5IO.hh
  HDF5Traits.hh
  HDF5Types.hh
  PyFileIO.hh
  SiloFileIO.hh
  vectorstringUtilities.hh
  )

set(FileIO_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/FileIO
  ${EIGEN_LIB_LOCATION}
  )


blt_add_library(NAME        Spheral_FileIO
                HEADERS     ${FileIO_headers}
                SOURCES     ${FileIO_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${FileIO_includes}
                )

install(TARGETS             Spheral_FileIO
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${FileIO_headers}
        DESTINATION include/FileIO
        )
list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_FileIO)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
