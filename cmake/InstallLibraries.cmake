# This function is used to force a build on a dependant project at cmake configuration phase.
# 
function (build_external_project target prefix url) #FOLLOWING ARGUMENTS are the CMAKE_ARGS of ExternalProject_Add
  set(trigger_build_dir ${CMAKE_BINARY_DIR}/force_${target})
  #mktemp dir in build tree
  file(MAKE_DIRECTORY ${trigger_build_dir} ${trigger_build_dir}/build)
  #generate false dependency project
  set(CMAKE_LIST_CONTENT "
    cmake_minimum_required(VERSION 3.10)

    include(ExternalProject)
    ExternalProject_add(${target}
      PREFIX ${prefix}/${target}
      GIT_REPOSITORY ${url}
      CMAKE_ARGS ${ARGN}
      )

    add_custom_target(trigger_${target})
    add_dependencies(trigger_${target} ${target})
  ")
  file(WRITE ${trigger_build_dir}/CMakeLists.txt "${CMAKE_LIST_CONTENT}")
  execute_process(COMMAND ${CMAKE_COMMAND} ..
    WORKING_DIRECTORY ${trigger_build_dir}/build
    )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . -j32
    WORKING_DIRECTORY ${trigger_build_dir}/build
    )
endfunction()


function (build_external_project_silo target prefix url) #FOLLOWING ARGUMENTS are the CMAKE_ARGS of ExternalProject_Add
  set(trigger_build_dir ${CMAKE_BINARY_DIR}/force_${target})
  #mktemp dir in build tree
  file(MAKE_DIRECTORY ${trigger_build_dir} ${trigger_build_dir}/build)
  #generate false dependency project
  set(CMAKE_LIST_CONTENT "
    cmake_minimum_required(VERSION 3.10)

    include(ExternalProject)
    ExternalProject_add(${target}
      PREFIX ${prefix}/${target}
      URL ${url}
      CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/../silo/silo/src/silo/configure --enable-shared=no --enable-fortran=no --with-hdf5=${PROJECT_SOURCE_DIR}/../hdf5/hdf5/include,${PROJECT_SOURCE_DIR}/../hdf5/hdf5/lib --prefix=${prefix}/${target} --enable-silex=no --enable-browser=yes
      BUILD_COMMAND make -j32
      INSTALL_COMMAND make install
      CMAKE_ARGS ${ARGN}
      )

    add_custom_target(trigger_${target})
    add_dependencies(trigger_${target} ${target})
  ")
  file(WRITE ${trigger_build_dir}/CMakeLists.txt "${CMAKE_LIST_CONTENT}")
  execute_process(COMMAND ${CMAKE_COMMAND} ..
    WORKING_DIRECTORY ${trigger_build_dir}/build
    )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . -j32
    WORKING_DIRECTORY ${trigger_build_dir}/build
    )
endfunction()


function (build_external_project_hdf5 target prefix url) #FOLLOWING ARGUMENTS are the CMAKE_ARGS of ExternalProject_Add
  set(trigger_build_dir ${CMAKE_BINARY_DIR}/force_${target})
  #mktemp dir in build tree
  file(MAKE_DIRECTORY ${trigger_build_dir} ${trigger_build_dir}/build)
  #generate false dependency project
  set(CMAKE_LIST_CONTENT "
    cmake_minimum_required(VERSION 3.10)

    include(ExternalProject)
    ExternalProject_add(${target}
      PREFIX ${prefix}/${target}
      URL ${url}
      CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/../hdf5/hdf5/src/hdf5/configure --prefix=${prefix}/${target}
      BUILD_COMMAND make -j32
      INSTALL_COMMAND make install
      CMAKE_ARGS ${ARGN}
      )

    add_custom_target(trigger_${target})
    add_dependencies(trigger_${target} ${target})
  ")
  file(WRITE ${trigger_build_dir}/CMakeLists.txt "${CMAKE_LIST_CONTENT}")
  execute_process(COMMAND ${CMAKE_COMMAND} ..
    WORKING_DIRECTORY ${trigger_build_dir}/build
    )
  execute_process(COMMAND ${CMAKE_COMMAND} --build . -j32
    WORKING_DIRECTORY ${trigger_build_dir}/build
    )
endfunction()


set(BOOST_MIN_VERSION "1.62.0")
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set(BOOST_DIR ${Boost_INCLUDE_DIRS}/..)
endif()


message("\n---------- BUILDING POLYTOPE ----------")
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/../polytope/include/polytope/polytope.hh)
  build_external_project(polytope 
    ${PROJECT_SOURCE_DIR}/../polytope/
    "https://github.com/pbtoast/polytope" 
    "-DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../polytope/"
    )
else()
  message(STATUS "Polytope Built")
endif()
set(POLYTOPE_DIR "${PROJECT_SOURCE_DIR}/../polytope")
message("---------------------------------------\n")


message("\n---------- BUILDING EIGEN ----------")
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/../eigen/include/eigen3/Eigen/Eigen)
  build_external_project(eigen 
    ${PROJECT_SOURCE_DIR}/../eigen/
    "https://github.com/eigenteam/eigen-git-mirror"
    "-DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../eigen/"
    )
else()
  message(STATUS "Eigen Built")
endif()
set(EIGEN_DIR "${PROJECT_SOURCE_DIR}/../eigen")
message("---------------------------------------\n")


message("\n---------- BUILDING QHULL ----------")
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/../qhull/lib/libqhullstatic.a)
  build_external_project(qhull
    ${PROJECT_SOURCE_DIR}/../qhull/
    "https://github.com/qhull/qhull"
    "-DCMAKE_C_FLAGS=-fPIC -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../qhull/"
    )
else()
  message(STATUS "Qhull Built")
endif()
set(QHULL_DIR "${PROJECT_SOURCE_DIR}/../qhull")
message("--------------------------------------\n")


message("\n---------- BUILDING PYBIND11 ----------")
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/../pybind11/include/pybind11/pybind11.h)
  build_external_project(pybind11
    ${PROJECT_SOURCE_DIR}/../pybind11/
    "https://github.com/pybind/pybind11"
    "-DPYBIND11_TEST=Off -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../pybind11/"
    )
else()
  message(STATUS "Pybind11 Built")
endif()
set(PYBIND11_DIR "${PROJECT_SOURCE_DIR}/../pybind11")
message("--------------------------------------\n")

message("\n---------- BUILDING HDF5 ----------")
#if(NOT EXISTS ${PROJECT_SOURCE_DIR}/../pybind11/include/pybind11/pybind11.h)
  build_external_project_hdf5(hdf5
    ${PROJECT_SOURCE_DIR}/../hdf5/
    "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.bz2
"
    "-DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../hdf5/"
    )
#else()
#  message(STATUS "Pybind11 Built")
#endif()
set(HDF5_DIR "${PROJECT_SOURCE_DIR}/../hdf5/hdf5")
message("--------------------------------------\n")

message("\n---------- BUILDING SILO ----------")
#if(NOT EXISTS ${PROJECT_SOURCE_DIR}/../pybind11/include/pybind11/pybind11.h)
  build_external_project_silo(silo
    ${PROJECT_SOURCE_DIR}/../silo/
    "https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz"
    "-DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../silo/"
    )
#else()
#  message(STATUS "Pybind11 Built")
#endif()
set(SILO_DIR "${PROJECT_SOURCE_DIR}/../silo/silo")
message("--------------------------------------\n")


