macro(spheral_add_executable)
  set(options )
  set(singleValueArgs NAME TEST REPRODUCER BENCHMARK)
  set(multiValueArgs SOURCES DEPENDS_ON)

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  if (ENABLE_OPENMP)
    list (APPEND arg_DEPENDS_ON openmp)
  endif ()

  if (ENABLE_CUDA)
    list (APPEND arg_DEPENDS_ON cuda)
  endif ()

  if (ENABLE_HIP)
    list (APPEND arg_DEPENDS_ON blt::hip)
    list (APPEND arg_DEPENDS_ON blt::hip_runtime)
  endif ()

  if (${arg_TEST})
    set (_output_dir ${CMAKE_BINARY_DIR}/test)
  elseif (${arg_REPRODUCER})
    set (_output_dir ${CMAKE_BINARY_DIR}/reproducers)
  elseif (${arg_BENCHMARK})
    set (_output_dir ${CMAKE_BINARY_DIR}/benchmark)
  else ()
    set (_output_dir ${CMAKE_BINARY_DIR}/bin)
  endif()

  blt_add_executable(
    NAME ${arg_NAME}
    SOURCES ${arg_SOURCES}
    DEPENDS_ON ${arg_DEPENDS_ON}
    OUTPUT_DIR ${_output_dir}
    )

  target_include_directories(${arg_NAME} SYSTEM PRIVATE ${SPHERAL_EXTERN_INCLUDES})
  target_include_directories(${arg_NAME} SYSTEM PRIVATE ${SPHERAL_ROOT_DIR}/src)
  target_include_directories(${arg_NAME} SYSTEM PRIVATE ${PROJECT_BINARY_DIR}/src)

endmacro(spheral_add_executable)

macro(spheral_add_test)
  set(options DEBUG_LINKER)
  set(singleValueArgs NAME)
  set(multiValueArgs SOURCES DEPENDS_ON)

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  set(original_test_name ${arg_NAME})
  set(original_src ${arg_SOURCES})
  set(original_deps ${arg_DEPENDS_ON})

  get_property(SPHERAL_BLT_DEPENDS GLOBAL PROPERTY SPHERAL_BLT_DEPENDS)

  blt_add_library(
    NAME ${original_test_name}_lib
    SOURCES ${TEST_LIB_SOURCE}
    SOURCES ${CMAKE_SOURCE_DIR}/src/spheralCXX.cc
    DEPENDS_ON ${SPHERAL_BLT_DEPENDS} ${original_deps}
    SHARED False 
    )

  target_link_options(${original_test_name}_lib PRIVATE "-Wl,--unresolved-symbols=ignore-in-object-files")

  spheral_add_executable(
    NAME ${original_test_name}.exe
    SOURCES ${original_src}
    DEPENDS_ON gtest ${CMAKE_THREAD_LIBS_INIT} ${original_test_name}_lib
    TEST On)

  blt_add_test(
    NAME ${original_test_name}
    COMMAND ${TEST_DRIVER} ${original_test_name}.exe)

  target_include_directories(${original_test_name}.exe SYSTEM PRIVATE ${SPHERAL_ROOT_DIR}/tests/cpp/include)

  if (${arg_DEBUG_LINKER})
    target_link_options(${original_test_name}.exe PUBLIC "-Wl,--warn-unresolved-symbols")
  else()
    target_link_options(${original_test_name}.exe PUBLIC "-Wl,--unresolved-symbols=ignore-all")
  endif()

endmacro(spheral_add_test)
