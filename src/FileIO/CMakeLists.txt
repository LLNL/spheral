set(FileIO_sources
    FileIO.cc
    FlatFileIO.cc
    PyFileIO.cc
    SiloFileIO.cc
    vectorstringUtilities.cc
    )

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

blt_add_library(NAME        Spheral_FileIO
                SOURCES     ${FileIO_sources}
                HEADERS     ${FileIO_headers}
                DEPENDS_ON  ${spheral_depends}
                )

install(TARGETS              Spheral_FileIO
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include
        )

install(FILES       ${FileIO_headers}
        DESTINATION include/FileIO
        )
