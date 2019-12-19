set(PATCH_DIR ${PROJECT_SOURCE_DIR}/src/thirdPartyLibs)

function(DownloadAndBuildLib TARGET_NAME)
  if(NOT EXISTS ${${TARGET_NAME}_EXISTS_FILE})
    set(${TARGET_NAME}_CMAKE_LIST_CONTENT "
      cmake_minimum_required(VERSION 3.10)
      include(ExternalProject)

      ${${TARGET_NAME}_EXTERNAL_PROJECT_FUNCTION}

      add_custom_target(trigger_${${TARGET_NAME}_TARGET})
      add_dependencies(trigger_${${TARGET_NAME}_TARGET} ${${TARGET_NAME}_TARGET})
    ")

    set(trigger_build_dir ${CMAKE_BINARY_DIR}/force_${TARGE_NAME})
    #mktemp dir in build tree
    file(MAKE_DIRECTORY ${trigger_build_dir} ${trigger_build_dir}/build)
    #generate false dependency project 
    file(WRITE ${trigger_build_dir}/CMakeLists.txt ${${TARGET_NAME}_CMAKE_LIST_CONTENT})
    execute_process(COMMAND ${CMAKE_COMMAND} ..
      WORKING_DIRECTORY ${trigger_build_dir}/build
    )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . -j32
      WORKING_DIRECTORY ${trigger_build_dir}/build
    )
    
    if(NOT EXISTS ${${TARGET_NAME}_EXISTS_FILE})
      message(ERROR "Failed to build ${TARGET_NAME}")
    endif()

  else()
    message(STATUS "${TARGET_NAME} Built")
  endif()
endfunction()



################################
# BOOST
################################
message("\n---------- BUILDING BOOST ----------")
set(BOOST_PREFIX ${PROJECT_SOURCE_DIR}/..)
set(BOOST_TARGET boost)
set(BOOST_DIR "${BOOST_PREFIX}/boost")
set(BOOST_EXISTS_FILE "${BOOST_PREFIX}/boost/lib/libboost_system.a")

#set(BOOST_URL "https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.bz2")
set(BOOST_URL "https://sourceforge.net/projects/boost/files/boost/1.69.0/boost_1_69_0.tar.bz2")
set(BOOST_SRC_DIR "${BOOST_PREFIX}/boost/src/boost")
set(BOOST_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${BOOST_TARGET}
    PREFIX ${BOOST_PREFIX}/${BOOST_TARGET}
    URL ${BOOST_URL} 
    #PATCH_COMMAND patch -t ${BOOST_SRC_DIR}/config/config.guess ${PATCH_DIR}/config.guess-boost-4.10.2-bsd.patch &&
    #              patch -t ${BOOST_SRC_DIR}/config/config.sub   ${PATCH_DIR}/config.sub-boost-4.10.2-bsd.patch
    CONFIGURE_COMMAND ${BOOST_SRC_DIR}/bootstrap.sh
    --with-toolset=${CMAKE_C_COMPILER_ID}
    --without-libraries=atomic,container,coroutine,log,chrono,context,date_time,exception,fiber,filesystem,graph,graph_parallel,iostreams,locale,math,mpi,program_options,python,random,regex,serialization,system,test,thread,timer,type_erasure,wave
    --prefix=${BOOST_DIR}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${BOOST_SRC_DIR}/b2 install
    INSTALL_COMMAND echo \"Skipping Install Step\"
  )
")
DownloadAndBuildLib(BOOST)
message("--------------------------------------\n")

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
# PYTHON
################################
message("\n---------- BUILDING PYTHON ----------")
set(PYTHON_PREFIX ${PROJECT_SOURCE_DIR}/..)
set(PYTHON_TARGET python)
set(PYTHON_DIR ${PYTHON_PREFIX}/python)
set(PYTHON_EXISTS_FILE "${PYTHON_DIR}/include/python2.7/Python.h")

set(PYTHON_URL "http://www.python.org/ftp/python/2.7.15/Python-2.7.15.tgz")
set(PYTHON_SRC_DIR "${PYTHON_PREFIX}/python/src/python")
set(PYTHON_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${PYTHON_TARGET}
    PREFIX ${PYTHON_PREFIX}/${PYTHON_TARGET}
    URL ${PYTHON_URL} 
    CONFIGURE_COMMAND ${PYTHON_SRC_DIR}/configure
                      --with-cxx-main='${CMAKE_CXX_COMPILER}'
                      --disable-ipv6
                      --prefix=${PYTHON_PREFIX}/${PYTHON_TARGET}
    BUILD_COMMAND make -j32
    INSTALL_COMMAND make install
  )
")

DownloadAndBuildLib(PYTHON)
message("--------------------------------------\n")



################################
# CONFUGURING PIP/PYTHON LIBS
################################
set(PYTHON_EXE ${PYTHON_DIR}/bin/python2.7)
set(PIP_EXE    ${PYTHON_DIR}/bin/pip2.7)
set(PYTHON_SITE_PACKAGE_DIR ${PYTHON_DIR}/lib/python2.7/site-packages)

# install pip and python packages
if (NOT EXISTS ${PIP_EXE})
  execute_process(COMMAND curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py)
  execute_process(COMMAND ${PYTHON_EXE} get-pip.py)
  execute_process(COMMAND ${PIP_EXE} install PYB11Generator -t ${PYTHON_SITE_PACKAGE_DIR})
  execute_process(COMMAND ${PIP_EXE} install mpi4py -t ${PYTHON_SITE_PACKAGE_DIR})
  execute_process(COMMAND ${PIP_EXE} install numpy -t ${PYTHON_SITE_PACKAGE_DIR})
endif()

# Find python
include(cmake/libraries/FindPython.cmake)



################################
# PYBIND11
################################
message("\n---------- BUILDING PYBIND11 ----------")
set(PYBIND11_PREFIX ${PROJECT_SOURCE_DIR}/../pybind11/)
set(PYBIND11_TARGET pybind11)
set(PYBIND11_DIR ${PYBIND11_PREFIX})
set(PYBIND11_EXISTS_FILE "${PYBIND11_PREFIX}/include/pybind11/pybind11.h")

set(PYBIND11_GIT "https://github.com/pybind/pybind11")
set(PYBIND11_CMAKE_ARGS "-DPYBIND11_TEST=Off -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../pybind11/")
set(PYBIND11_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${PYBIND11_TARGET}
    PREFIX ${PYBIND11_PREFIX}/${PYBIND11_TARGET}
    GIT_REPOSITORY ${PYBIND11_GIT}
    CMAKE_ARGS ${PYBIND11_CMAKE_ARGS}
  )
")
DownloadAndBuildLib(PYBIND11)
message("--------------------------------------\n")

if (PYBIND11_DIR)
  include(cmake/libraries/FindPybind11.cmake)
  if (PYBIND11_FOUND)
    blt_register_library( NAME      pybind11
                          INCLUDES  ${PYBIND11_INCLUDE_DIRS}
                          #LIBRARIES ${PYBIND11_LIBRARY}
                          DEFINES   HAVE_PYBIND11=1 )
  else()
    message(FATAL_ERROR "Unable to find PYBIND11 with given path. ${PYBIND11_DIR}")
  endif()
else()
  message(STATUS "Library Disabled: PYBIND11")
  set(BLT_PYBIND11_DEFINES "HAVE_PYBIND11=0" CACHE PATH "")
endif()



################################
# POLYTOPE
################################
message("\n---------- BUILDING POLYTOPE ----------")
set(POLYTOPE_PREFIX ${PROJECT_SOURCE_DIR}/../polytope/)
set(POLYTOPE_TARGET polytope)
set(POLYTOPE_DIR ${POLYTOPE_PREFIX})
set(POLYTOPE_EXISTS_FILE "${POLYTOPE_PREFIX}/include/polytope/polytope.hh")

set(POLYTOPE_GIT "https://github.com/pbtoast/polytope")
set(POLYTOPE_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${POLYTOPE_PREFIX} 
                         -DPYBIND11_INCLUDE_DIRS=${PYBIND11_INCLUDE_DIRS}
                         -DUSE_PYTHON=On
                         -DPYTHON_EXE=${PYTHON_EXE}
                         -DBoost_INCLUDE_DIR=${BOOST_INCLUDE_DIRS}
                         -DTESTING=Off")
set(POLYTOPE_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${POLYTOPE_TARGET}
    PREFIX ${POLYTOPE_PREFIX}/${POLYTOPE_TARGET}
    PATCH_COMMAND patch -t ${POLYTOPE_PREFIX}/polytope/src/polytope/CMakeLists.txt            ${PATCH_DIR}/polytope-CMakeLists.patch &&
                  patch -t ${POLYTOPE_PREFIX}/polytope/src/polytope/src/PYB11/CMakeLists.txt  ${PATCH_DIR}/polytope-PYB11-CMakeLists.patch

    GIT_REPOSITORY ${POLYTOPE_GIT}
    CMAKE_ARGS ${POLYTOPE_CMAKE_ARGS}
  )
")
DownloadAndBuildLib(POLYTOPE)
message("---------------------------------------\n")

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



################################
# EIGEN
################################
message("\n---------- BUILDING EIGEN ----------")
set(EIGEN_PREFIX ${PROJECT_SOURCE_DIR}/../eigen/)
set(EIGEN_TARGET eigen)
set(EIGEN_DIR ${EIGEN_PREFIX})
set(EIGEN_EXISTS_FILE "${EIGEN_PREFIX}/include/eigen3/Eigen/Eigen")

set(EIGEN_GIT "https://github.com/eigenteam/eigen-git-mirror")
set(EIGEN_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${EIGEN_PREFIX}")
set(EIGEN_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${EIGEN_TARGET}
    PREFIX ${EIGEN_PREFIX}/${EIGEN_TARGET}
    GIT_REPOSITORY ${EIGEN_GIT}
    CMAKE_ARGS ${EIGEN_CMAKE_ARGS}
  )
")
DownloadAndBuildLib(EIGEN)
message("---------------------------------------\n")

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
message("\n---------- BUILDING QHULL ----------")
set(QHULL_PREFIX ${PROJECT_SOURCE_DIR}/..)
set(QHULL_TARGET qhull)
set(QHULL_DIR ${QHULL_PREFIX}/qhull)
set(QHULL_EXISTS_FILE "${QHULL_DIR}/lib/libqhullstatic.a")

set(QHULL_GIT "https://github.com/qhull/qhull")
#set(QHULL_CMAKE_ARGS "-DCMAKE_C_FLAGS=-fPIC -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../qhull/")
set(QHULL_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_C_FLAGS=-fPIC -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${QHULL_DIR}")
set(QHULL_SRC_DIR ${QHULL_DIR}/src/qhull/src)
set(QHULL_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${QHULL_TARGET}
    PREFIX ${QHULL_PREFIX}/${QHULL_TARGET}
    PATCH_COMMAND patch -t ${QHULL_SRC_DIR}/libqhull/qhull_a.h ${PATCH_DIR}/qhull-2015.2-qhull_a.h-patch &&
                  patch -t ${QHULL_SRC_DIR}/libqhull_r/qhull_ra.h ${PATCH_DIR}/qhull-2015.2-qhull_ra.h-patch
    GIT_REPOSITORY ${QHULL_GIT}
    CMAKE_ARGS ${QHULL_CMAKE_ARGS}
  )
")
DownloadAndBuildLib(QHULL)
message("--------------------------------------\n")

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
# HDF5
################################
message("\n---------- BUILDING HDF5 ----------")
set(HDF5_PREFIX ${PROJECT_SOURCE_DIR}/..)
set(HDF5_TARGET hdf5)
set(HDF5_DIR "${HDF5_PREFIX}/hdf5")
set(HDF5_EXISTS_FILE "${HDF5_PREFIX}/hdf5/lib/libhdf5.a")

set(HDF5_URL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.bz2")
set(HDF5_SRC_DIR "${HDF5_PREFIX}/hdf5/src/hdf5")
set(HDF5_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${HDF5_TARGET}
    PREFIX ${HDF5_PREFIX}/${HDF5_TARGET}
    URL ${HDF5_URL} 
    CONFIGURE_COMMAND ${HDF5_SRC_DIR}/configure 
                      --prefix=${HDF5_PREFIX}/${HDF5_TARGET}
    BUILD_COMMAND make -j32
    INSTALL_COMMAND make install
  )
")
DownloadAndBuildLib(HDF5)
message("--------------------------------------\n")

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
# SILO
################################
message("\n---------- BUILDING SILO ----------")
set(SILO_PREFIX ${PROJECT_SOURCE_DIR}/..)
set(SILO_TARGET silo)
set(SILO_DIR "${SILO_PREFIX}/silo")
set(SILO_EXISTS_FILE "${SILO_PREFIX}/silo/lib/libsiloh5.a")

set(SILO_URL "https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz")
set(SILO_SRC_DIR "${SILO_PREFIX}/silo/src/silo")
set(SILO_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${SILO_TARGET}
    PREFIX ${SILO_PREFIX}/${SILO_TARGET}
    URL ${SILO_URL} 
    PATCH_COMMAND patch -t ${SILO_SRC_DIR}/config/config.guess ${PATCH_DIR}/config.guess-silo-4.10.2-bsd.patch &&
                  patch -t ${SILO_SRC_DIR}/config/config.sub   ${PATCH_DIR}/config.sub-silo-4.10.2-bsd.patch
    CONFIGURE_COMMAND ${SILO_SRC_DIR}/configure
                      --enable-shared=no
                      --enable-fortran=no
                      --with-hdf5=${HDF5_DIR}/include,${HDF5_DIR}/lib
                      --prefix=${SILO_PREFIX}/${SILO_TARGET}
                      --enable-silex=no
                      --enable-browser=yes
    BUILD_COMMAND make -j32
    INSTALL_COMMAND make install
  )
")
DownloadAndBuildLib(SILO)
message("--------------------------------------\n")

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
