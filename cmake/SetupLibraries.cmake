################################
# SILO
################################
if (SILO_DIR)
    include(cmake/libraries/FindSILO.cmake)
    if (SILO_FOUND)
        blt_register_library( NAME       silo
                              INCLUDES   ${SILO_INCLUDE_DIRS}
                              LIBRARIES  ${SILO_LIBRARY}
                              DEPENDS_ON hdf5
                              DEFINES    HAVE_SILO=1 )
    else()
        message(FATAL_ERROR "Unable to find SILO with given path.")
    endif()
else()
    message(STATUS "Library Disabled: SILO")
    set(BLT_SILO_DEFINES "HAVE_SILO=0" CACHE PATH "")
endif()

################################
# HDF5
################################
if (HDF5_DIR)
    include(cmake/libraries/FindHDF5.cmake)
    if (HDF5_FOUND)
        blt_register_library( NAME       hdf5
                              INCLUDES   ${HDF5_INCLUDE_DIRS}
                              LIBRARIES  ${HDF5_LIBRARIES}
                              DEFINES    HAVE_HDF5=1 )
    else()
        message(FATAL_ERROR "Unable to find HDF5 with given path.")
    endif()
else()
    message(STATUS "Library Disabled: HDF5")
    set(BLT_HDF5_DEFINES "HAVE_HDF5=0" CACHE PATH "")
endif()

################################
# BOOST
################################
if (BOOST_DIR)
    include(cmake/libraries/FindBOOST.cmake)
    if (BOOST_FOUND)
        blt_register_library( NAME       BOOST
                              INCLUDES   ${BOOST_INCLUDE_DIRS}
                              DEFINES    HAVE_BOOST=1 )
    else()
        message(FATAL_ERROR "Unable to find BOOST with given path.")
    endif()
else()
    message(STATUS "Library Disabled: BOOST")
    set(BLT_BOOST_DEFINES "HAVE_BOOST=0" CACHE PATH "")
endif()

################################
# EIGEN
################################
if (EIGEN_DIR)
    include(cmake/libraries/FindEigen.cmake)
    if (EIGEN_FOUND)
        blt_register_library( NAME       eigen
                              INCLUDES   ${EIGEN_INCLUDE_DIRS}
                              DEFINES    HAVE_EIGEN=1 )
    else()
        message(FATAL_ERROR "Unable to find EIGEN with given path.")
    endif()
else()
    message(STATUS "Library Disabled: EIGEN")
    set(BLT_EIGEN_DEFINES "HAVE_EIGEN=0" CACHE PATH "")
endif()

################################
# QHULL
################################
if (QHULL_DIR)
    include(cmake/libraries/FindQhull.cmake)
    if (QHULL_FOUND)
        blt_register_library( NAME       qhull
                              INCLUDES   ${QHULL_INCLUDE_DIRS}
                              LIBRARIES  ${QHULL_LIBRARY}
                              DEFINES    HAVE_QHULL=1 )
    else()
        message(FATAL_ERROR "Unable to find QHULL with given path.")
    endif()
else()
    message(STATUS "Library Disabled: QHULL")
    set(BLT_QHULL_DEFINES "HAVE_QHULL=0" CACHE PATH "")
endif()

################################
# POLYTOPE
################################
if (POLYTOPE_DIR)
    include(cmake/libraries/FindPolytope.cmake)
    if (POLYTOPE_FOUND)
        blt_register_library( NAME       polytope
                              INCLUDES   ${POLYTOPE_INCLUDE_DIRS}
                              LIBRARIES  ${POLYTOPE_LIBRARY}
                              DEFINES    HAVE_POLYTOPE=1 )
    else()
        message(FATAL_ERROR "Unable to find POLYTOPE with given path.")
    endif()
else()
    message(STATUS "Library Disabled: POLYTOPE")
    set(BLT_POLYTOPE_DEFINES "HAVE_POLYTOPE=0" CACHE PATH "")
endif()
