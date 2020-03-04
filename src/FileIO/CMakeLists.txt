include_directories(.)
set(FileIO_sources
    FileIO.cc
    FlatFileIO.cc
    SiloFileIO.cc
    PyFileIO.cc
    vectorstringUtilities.cc
    )

set(FileIO_headers
    DbFileIO.hh
    FileIO.hh
    FileIOInline.hh
    HDF5IO.hh
    HDF5Traits.hh
    HDF5Types.hh
    )

set(FileIO_python_targets
  GzipFileIO.py
  PlyFileIO.py
  )

spheral_add_cxx_library(FileIO)

install(TARGETS              Spheral_FileIO
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${FileIO_headers}
        DESTINATION include/FileIO
        )
      
if(NOT ENABLE_CXXONLY)
  install(FILES       ${FileIO_python_targets}
          DESTINATION .
          )
endif()

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_FileIO)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
