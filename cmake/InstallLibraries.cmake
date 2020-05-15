set(PATCH_DIR ${PROJECT_SOURCE_DIR}/thirdPartyLibs)
set(CACHE_DIR ${PROJECT_SOURCE_DIR}/thirdPartyLibs/cache)

if(SPHERAL_TPL_DIR)
  if(NOT EXISTS SPHERAL_TPL_DIR)
    get_filename_component(SPHERAL_TPL_DIR_PARENT ${SPHERAL_TPL_DIR} DIRECTORY)
    if(NOT EXISTS "${SPHERAL_TPL_DIR_PARENT}/")
      set(SPHERAL_TPL_DIR ${PROJECT_SOURCE_DIR}/..)
      message(FATAL_ERROR "No Directory ${SPHERAL_TPL_DIR_PARENT}/")
    endif()
  endif()
else()
  set(SPHERAL_TPL_DIR ${PROJECT_SOURCE_DIR}/..)
endif()
set(SPHERAL_TPL_DIR ${SPHERAL_TPL_DIR}/tpl)
message(STATUS "TPL Install Dir at: ${SPHERAL_TPL_DIR}")



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
  execute_process(COMMAND ${CMAKE_COMMAND} --build . 
      WORKING_DIRECTORY ${trigger_build_dir}/build
    )
    
    if(NOT EXISTS ${${TARGET_NAME}_EXISTS_FILE})
      message(FATAL_ERROR "Failed to build ${TARGET_NAME}")
    endif()

  else()
    message(STATUS "${TARGET_NAME} Built")
  endif()
endfunction()



################################
# BOOST
################################
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(TOOLSET "gcc")
else()
  set(TOOLSET ${CMAKE_CXX_COMPILER_ID})
endif()

if(INSTALL_TPLS AND NOT BOOST_DIR)
  message("\n---------- BUILDING BOOST ----------")
  set(BOOST_PREFIX ${SPHERAL_TPL_DIR})
  set(BOOST_TARGET boost)
  set(BOOST_DIR "${BOOST_PREFIX}/boost")
  set(BOOST_EXISTS_FILE "${BOOST_PREFIX}/boost/lib/libboost_system.a")
  set(BOOST_CACHE "${CACHE_DIR}/boost_1_69_0.tar.bz2")

  if (NOT BOOST_URL)
    if (EXISTS ${BOOST_CACHE})
      set(BOOST_URL ${BOOST_CACHE})
    else()
      set(BOOST_URL "https://sourceforge.net/projects/boost/files/boost/1.69.0/boost_1_69_0.tar.bz2")
    endif()

  endif()

  set(BOOST_SRC_DIR "${BOOST_PREFIX}/boost/src/boost")
  set(BOOST_EXTERNAL_PROJECT_FUNCTION "
    ExternalProject_add(${BOOST_TARGET}
      PREFIX ${BOOST_PREFIX}/${BOOST_TARGET}
      URL ${BOOST_URL} 
      DOWNLOAD_DIR ${CACHE_DIR}
      #PATCH_COMMAND patch -t ${BOOST_SRC_DIR}/config/config.guess ${PATCH_DIR}/config.guess-boost-4.10.2-bsd.patch &&
      #              patch -t ${BOOST_SRC_DIR}/config/config.sub   ${PATCH_DIR}/config.sub-boost-4.10.2-bsd.patch
      CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${BOOST_SRC_DIR}/bootstrap.sh
      --with-toolset=${TOOLSET}
      --without-libraries=atomic,container,coroutine,log,chrono,context,date_time,exception,fiber,filesystem,graph,graph_parallel,iostreams,locale,math,mpi,program_options,python,random,regex,serialization,system,test,thread,timer,type_erasure,wave
      --prefix=${BOOST_DIR}
      BUILD_IN_SOURCE 1
      BUILD_COMMAND ${BOOST_SRC_DIR}/b2 install
      INSTALL_COMMAND echo \"Skipping Install Step\"
    )
  ")
  DownloadAndBuildLib(BOOST)
  message("--------------------------------------\n")
endif()

if (BOOST_DIR)
    include(../cmake/libraries/FindBOOST.cmake)
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
if(INSTALL_TPLS AND NOT PYTHON_DIR)
  message("\n---------- BUILDING PYTHON ----------")
  get_filename_component(PYTHON_PREFIX ${SPHERAL_TPL_DIR} DIRECTORY)
  set(PYTHON_PREFIX "${PYTHON_PREFIX}")
  set(PYTHON_TARGET python)
  set(PYTHON_DIR ${PYTHON_PREFIX}/python)
  set(PYTHON_EXISTS_FILE "${PYTHON_DIR}/include/python2.7/Python.h")
  set(PYTHON_CACHE "${CACHE_DIR}/Python-2.7.15.tgz")

  if(NOT PYTHON_URL)
    if (EXISTS ${PYTHON_CACHE})
      set(PYTHON_URL ${PYTHON_CACHE})
    else()
      set(PYTHON_URL "http://www.python.org/ftp/python/2.7.15/Python-2.7.15.tgz")
    endif()
  endif()

  set(_PYTHON_SRC_DIR "${PYTHON_PREFIX}/python/src/python")
  set(PYTHON_EXTERNAL_PROJECT_FUNCTION "
    ExternalProject_add(${PYTHON_TARGET}
      PREFIX ${PYTHON_PREFIX}/${PYTHON_TARGET}
      URL ${PYTHON_URL} 
      DOWNLOAD_DIR ${CACHE_DIR}
      CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${_PYTHON_SRC_DIR}/configure
                        --with-cxx-main='${CMAKE_CXX_COMPILER}'
                        --disable-ipv6
                        --prefix=${PYTHON_PREFIX}/${PYTHON_TARGET}
      BUILD_COMMAND make 
      INSTALL_COMMAND make install
    )
  ")

  DownloadAndBuildLib(PYTHON)
  message("--------------------------------------\n")
endif()

set(PIPTARGETS 
  setuptools
  PYB11Generator
  mpi4py
  numpy
  numpy-stl
  matplotlib
  decorator
  h5py
  sphinx
  sphinx_rtd_theme
  twine
  cython
  sobol
  scipy
  pipreqs
  )

if(PYTHON_DIR)
  ################################
  # CONFUGURING PIP/PYTHON LIBS
  ################################
  set(PYTHON_EXE ${PYTHON_DIR}/bin/python2.7)
  set(PYTHON_SITE_PACKAGE_DIR ${PYTHON_DIR}/lib/python2.7/site-packages)
  set(PYTHON_EXECUTABLE ${PYTHON_EXE})
  set(PIP_EXE    ${PYTHON_DIR}/bin/pip2.7)
  set(PYTHON_VERSION "2.7")

  if(NOT ENABLE_CXXONLY)
    include(FindPythonModule)

    find_python_module(pip QUIET)
    if(NOT PY_PIP)
      set(PIP_DIST pip-9.0.1-py2.py3-none-any.whl)
      if (NOT EXISTS ${CACHE_DIR}/${PIP_DIST})
        execute_process(COMMAND wget https://pypi.python.org/packages/b6/ac/7015eb97dc749283ffdec1c3a88ddb8ae03b8fad0f0e611408f196358da3/pip-9.0.1-py2.py3-none-any.whl#md5=297dbd16ef53bcef0447d245815f5144
  -O ${CACHE_DIR}/${PIP_DIST})
    endif()
      execute_process(COMMAND ${PYTHON_EXE} ${CACHE_DIR}/${PIP_DIST}/pip install --no-index ${CACHE_DIR}/${PIP_DIST})
    endif()

    # Download and install PIP modules
    foreach(_target ${PIPTARGETS})
      # Case for numpy-stl
      set(_tmp_target ${_target})
      if(${_target} STREQUAL "numpy-stl")
        set(_tmp_target stl)
      endif()

      string(TOUPPER ${_tmp_target} _target_upper)
      find_python_module(${_tmp_target} QUIET)

      if (NOT PY_${_target_upper})
        execute_process(COMMAND ${PIP_EXE} download --no-binary :all -d ${CACHE_DIR} ${_target})
      endif()

      if (NOT PY_${_target_upper})
        execute_process(COMMAND ${PIP_EXE} install --upgrade ${_target} --no-index --find-links ${CACHE_DIR})
      endif()
    endforeach()

    find_python_module(Gnuplot QUIET)
    if(NOT PY_GNUPLOT)
      execute_process(COMMAND wget http://downloads.sourceforge.net/gnuplot-py/gnuplot-py-1.8.tar.gz -O ${CACHE_DIR}/gnuplot-py-1.8.tar.gz)
      execute_process(COMMAND tar -xvf ${CACHE_DIR}/gnuplot-py-1.8.tar.gz)
      execute_process(COMMAND ${PYTHON_EXE} setup.py install
                      WORKING_DIRECTORY gnuplot-py-1.8)
    endif()
  endif()
endif()

# Find python
include(../cmake/libraries/FindPython.cmake)


if (NOT ENABLE_CXXONLY)
################################
# PYBIND11
################################
if(INSTALL_TPLS AND NOT PYBIND11_DIR AND NOT ENABLE_CXXONLY)
  message("\n---------- BUILDING PYBIND11 ----------")
  set(PYBIND11_PREFIX ${SPHERAL_TPL_DIR}/pybind11/)
  set(PYBIND11_TARGET pybind11)
  set(PYBIND11_DIR ${PYBIND11_PREFIX})
  set(PYBIND11_EXISTS_FILE "${PYBIND11_PREFIX}/include/pybind11/pybind11.h")
  set(PYBIND11_CACHE "${CACHE_DIR}/v2.4.3.tar.gz")

  if (NOT PYBIND11_URL)
    if (EXISTS ${PYBIND11_CACHE})
      set(PYBIND11_URL ${PYBIND11_CACHE})
    else()
      set(PYBIND11_URL "https://github.com/pybind/pybind11/archive/v2.4.3.tar.gz")
    endif()
  endif()

  set(PYBIND11_CMAKE_ARGS "-DPYBIND11_TEST=Off -DCMAKE_INSTALL_PREFIX=${PYBIND11_PREFIX} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} ")
  set(PYBIND11_EXTERNAL_PROJECT_FUNCTION "
    ExternalProject_add(${PYBIND11_TARGET}
      PREFIX ${PYBIND11_PREFIX}/${PYBIND11_TARGET}
      URL ${PYBIND11_URL}
      DOWNLOAD_DIR ${CACHE_DIR}
      CMAKE_ARGS ${PYBIND11_CMAKE_ARGS}
    )
  ")
  DownloadAndBuildLib(PYBIND11)
  message("--------------------------------------\n")
endif()

if (PYBIND11_DIR)
  include(../cmake/libraries/FindPybind11.cmake)
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

endif()


################################
# POLYTOPE
################################
if(INSTALL_TPLS AND NOT POLYTOPE_DIR)
  message("\n---------- BUILDING POLYTOPE ----------")
  set(POLYTOPE_PREFIX ${SPHERAL_TPL_DIR}/polytope/)
  set(POLYTOPE_TARGET polytope)
  set(POLYTOPE_DIR ${POLYTOPE_PREFIX})
  set(POLYTOPE_EXISTS_FILE "${POLYTOPE_PREFIX}/include/polytope/polytope.hh")
  set(POLYTOPE_CACHE "${CACHE_DIR}/0.6.2.tar.gz")

  if (NOT POLYTOPE_URL)
    if (EXISTS ${POLYTOPE_CACHE})
      set(POLYTOPE_URL ${POLYTOPE_CACHE})
    else()
      set(POLYTOPE_URL "https://github.com/pbtoast/polytope/archive/0.6.2.tar.gz")
    endif()
  endif()

  set(POLYTOPE_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                           -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                           -DCMAKE_INSTALL_PREFIX=${POLYTOPE_PREFIX} 
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

      URL ${POLYTOPE_URL}
      DOWNLOAD_DIR ${CACHE_DIR}
      CMAKE_ARGS ${POLYTOPE_CMAKE_ARGS}
    )
  ")
  DownloadAndBuildLib(POLYTOPE)
  message("---------------------------------------\n")
endif()

if (POLYTOPE_DIR)
    include(../cmake/libraries/FindPolytope.cmake)
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
if(INSTALL_TPLS AND NOT EIGEN_DIR)
  message("\n---------- BUILDING EIGEN ----------")
  set(EIGEN_PREFIX ${SPHERAL_TPL_DIR}/eigen/)
  set(EIGEN_TARGET eigen)
  set(EIGEN_DIR ${EIGEN_PREFIX})
  set(EIGEN_EXISTS_FILE "${EIGEN_PREFIX}/include/eigen3/Eigen/Eigen")
  set(EIGEN_CACHE "${CACHE_DIR}/3.3.7.tar.gz")

  if (NOT EIGEN_URL)
    if (EXISTS ${EIGEN_CACHE})
      set(EIGEN_URL ${EIGEN_CACHE})
    else()
      set(EIGEN_URL "https://github.com/eigenteam/eigen-git-mirror/archive/3.3.7.tar.gz")
    endif()
  endif()

  set(EIGEN_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_INSTALL_PREFIX=${EIGEN_PREFIX}")
  set(EIGEN_EXTERNAL_PROJECT_FUNCTION "
    ExternalProject_add(${EIGEN_TARGET}
      PREFIX ${EIGEN_PREFIX}/${EIGEN_TARGET}
      URL ${EIGEN_URL}
      DOWNLOAD_DIR ${CACHE_DIR}
      CMAKE_ARGS ${EIGEN_CMAKE_ARGS}
    )
  ")
  DownloadAndBuildLib(EIGEN)
  message("---------------------------------------\n")
endif()

if (EIGEN_DIR)
    include(../cmake/libraries/FindEigen.cmake)
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
if(INSTALL_TPLS AND NOT QHULL_DIR)
  message("\n---------- BUILDING QHULL ----------")
  set(QHULL_PREFIX ${SPHERAL_TPL_DIR})
  set(QHULL_TARGET qhull)
  set(QHULL_DIR ${QHULL_PREFIX}/qhull)
  set(QHULL_EXISTS_FILE "${QHULL_DIR}/lib/libqhullstatic.a")
  set(QHULL_CACHE "${CACHE_DIR}/2019.1.tar.gz")

  if (NOT QHULL_URL)
    if (EXISTS ${QHULL_CACHE})
      set(QHULL_URL ${QHULL_CACHE})
    else()
      set(QHULL_URL "https://github.com/qhull/qhull/archive/2019.1.tar.gz")
    endif()
  endif()

  set(QHULL_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_FLAGS=-fPIC -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${QHULL_DIR}")
  set(QHULL_SRC_DIR ${QHULL_DIR}/src/qhull/src)
  set(QHULL_EXTERNAL_PROJECT_FUNCTION "
    ExternalProject_add(${QHULL_TARGET}
      PREFIX ${QHULL_PREFIX}/${QHULL_TARGET}
      PATCH_COMMAND patch -t ${QHULL_SRC_DIR}/libqhull/qhull_a.h ${PATCH_DIR}/qhull-2015.2-qhull_a.h-patch &&
                    patch -t ${QHULL_SRC_DIR}/libqhull_r/qhull_ra.h ${PATCH_DIR}/qhull-2015.2-qhull_ra.h-patch
      URL ${QHULL_URL}
      DOWNLOAD_DIR ${CACHE_DIR}
      CMAKE_ARGS ${QHULL_CMAKE_ARGS}
  )
")
  DownloadAndBuildLib(QHULL)
  message("--------------------------------------\n")
endif()

if (QHULL_DIR)
    include(../cmake/libraries/FindQhull.cmake)
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
if(INSTALL_TPLS AND NOT HDF5_DIR)
  message("\n---------- BUILDING HDF5 ----------")
  set(HDF5_PREFIX ${SPHERAL_TPL_DIR})
  set(HDF5_TARGET hdf5)
  set(HDF5_DIR "${HDF5_PREFIX}/hdf5")
  set(HDF5_EXISTS_FILE "${HDF5_PREFIX}/hdf5/lib/libhdf5.a")
  set(HDF5_CACHE "${CACHE_DIR}/hdf5-1.10.4.tar.bz2")

  if (NOT HDF5_URL)
    if (EXISTS ${HDF5_CACHE})
      set(HDF5_URL ${HDF5_CACHE})
    else()
      set(HDF5_URL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.bz2")
    endif()
  endif()

  set(HDF5_SRC_DIR "${HDF5_PREFIX}/hdf5/src/hdf5")
  set(HDF5_EXTERNAL_PROJECT_FUNCTION "
    ExternalProject_add(${HDF5_TARGET}
      PREFIX ${HDF5_PREFIX}/${HDF5_TARGET}
      URL ${HDF5_URL} 
      DOWNLOAD_DIR ${CACHE_DIR}
      CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${HDF5_SRC_DIR}/configure 
                        --prefix=${HDF5_PREFIX}/${HDF5_TARGET}
      BUILD_COMMAND make 
      INSTALL_COMMAND make install
    )
  ")
  DownloadAndBuildLib(HDF5)
  message("--------------------------------------\n")
endif()

if (HDF5_DIR)
    include(../cmake/libraries/FindHDF5.cmake)
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
if(INSTALL_TPLS AND NOT SILO_DIR)
  message("\n---------- BUILDING SILO ----------")
  set(SILO_PREFIX ${SPHERAL_TPL_DIR})
  set(SILO_TARGET silo)
  set(SILO_DIR "${SILO_PREFIX}/silo")
  set(SILO_EXISTS_FILE "${SILO_PREFIX}/silo/lib/libsiloh5.a")
  set(SILO_CACHE "${CACHE_DIR}/silo-4.10.2-bsd.tar.gz")

  if(NOT SILO_URL)
    if (EXISTS ${SILO_CACHE})
      set(SILO_URL ${SILO_CACHE})
    else()
      set(SILO_URL "https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz")
    endif()
  endif()

  set(SILO_SRC_DIR "${SILO_PREFIX}/silo/src/silo")
  set(SILO_EXTERNAL_PROJECT_FUNCTION "
    ExternalProject_add(${SILO_TARGET}
      PREFIX ${SILO_PREFIX}/${SILO_TARGET}
      URL ${SILO_URL} 
      DOWNLOAD_DIR ${CACHE_DIR}
      PATCH_COMMAND patch -t ${SILO_SRC_DIR}/config/config.guess ${PATCH_DIR}/config.guess-silo-4.10.2-bsd.patch &&
                    patch -t ${SILO_SRC_DIR}/config/config.sub   ${PATCH_DIR}/config.sub-silo-4.10.2-bsd.patch
      CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${SILO_SRC_DIR}/configure
                        --enable-shared=no
                        --enable-fortran=no
                        --with-hdf5=${HDF5_DIR}/include,${HDF5_DIR}/lib
                        --prefix=${SILO_PREFIX}/${SILO_TARGET}
                        --enable-silex=no
                        --enable-browser=yes
      BUILD_COMMAND make 
      INSTALL_COMMAND make install
    )
  ")
  DownloadAndBuildLib(SILO)
  message("--------------------------------------\n")
endif()

if (SILO_DIR)
    include(../cmake/libraries/FindSILO.cmake)
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
