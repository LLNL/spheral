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
#if(NOT EXISTS ${PROJECT_SOURCE_DIR}/../qhull/lib/libqhullstatic.a)
  build_external_project(pybind11
    ${PROJECT_SOURCE_DIR}/../pybind11/
    "https://github.com/pybind/pybind11"
    "-DPYBIND11_TEST=Off -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/../pybind11/"
    )
  #else()
  #  message(STATUS "Qhull Built")
  #endif()
set(PYBIND11_DIR "${PROJECT_SOURCE_DIR}/../pybind11")
message("--------------------------------------\n")
