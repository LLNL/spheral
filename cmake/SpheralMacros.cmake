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

endmacro(spheral_add_executable)

macro(spheral_add_test)
  set(options )
  set(singleValueArgs NAME)
  set(multiValueArgs SOURCES DEPENDS_ON)

  cmake_parse_arguments(arg
    "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

  list (APPEND arg_DEPENDS_ON gtest Spheral_CXX ${CMAKE_THREAD_LIBS_INIT})

  set(original_test_name ${arg_NAME})

  spheral_add_executable(
    NAME ${arg_NAME}.exe
    SOURCES ${arg_SOURCES}
    DEPENDS_ON ${arg_DEPENDS_ON}
    TEST On)

  blt_add_test(
    NAME ${arg_NAME}
    #COMMAND ${TEST_DRIVER} $<TARGET_FILE:${arg_NAME}>)
    COMMAND ${TEST_DRIVER} ${arg_NAME})

  #spheral_set_failtest(${original_test_name})
endmacro(spheral_add_test)

