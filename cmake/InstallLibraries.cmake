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



set(BOOST_MIN_VERSION "1.62.0")
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set(BOOST_DIR ${Boost_INCLUDE_DIRS}/..)
endif()


message("\n---------- BUILDING POLYTOPE ----------")
set(POLYTOPE_PREFIX ${PROJECT_SOURCE_DIR}/../polytope/)
set(POLYTOPE_TARGET polytope)
set(POLYTOPE_DIR ${POLYTOPE_PREFIX})
set(POLYTOPE_EXISTS_FILE "${POLYTOPE_PREFIX}/include/polytope/polytope.hh")

set(POLYTOPE_GIT "https://github.com/pbtoast/polytope")
set(POLYTOPE_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${POLYTOPE_PREFIX}")
set(POLYTOPE_EXTERNAL_PROJECT_FUNCTION "
  ExternalProject_add(${POLYTOPE_TARGET}
    PREFIX ${POLYTOPE_PREFIX}/${POLYTOPE_TARGET}
    GIT_REPOSITORY ${POLYTOPE_GIT}
    CMAKE_ARGS ${POLYTOPE_CMAKE_ARGS}
  )
")
DownloadAndBuildLib(POLYTOPE)
message("---------------------------------------\n")


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
