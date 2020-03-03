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

#include(ExternalProject)
#ExternalProject_Add(polytope
#  GIT_REPOSITORY https://github.com/pbtoast/polytope
#  CMAKE_ARGS     "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/tpl-install"
#  )

message("\n---------- BUILDING POLYTOPE ----------")
if(NOT EXISTS ${PROJECT_BINARY_DIR}/tpl-install/include/polytope/polytope.hh)
  build_external_project(polytope 
    ${PROJECT_BINARY_DIR}/tpl-install/ 
    "https://github.com/pbtoast/polytope" 
    "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/tpl-install"
    )
else()
  message(STATUS "Polytope Built")
endif()
list(APPEND spheral_depends polytope)
message("---------------------------------------\n")

message("\n---------- BUILDING QHULL ----------")
if(NOT EXISTS ${PROJECT_BINARY_DIR}/tpl-install/lib/libqhullstatic.a)
  build_external_project(qhull
    ${PROJECT_BINARY_DIR}/tpl-install/
    "https://github.com/qhull/qhull"
    "-DCMAKE_C_FLAGS=-fPIC -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/tpl-install"
    )
else()
  message(STATUS "Qhull Built")
endif()
list(APPEND spheral_depends polytope)
message("--------------------------------------\n")



set(spheral_depends "${spheral_depends}" PARENT_SCOPE)
