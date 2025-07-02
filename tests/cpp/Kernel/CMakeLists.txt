include(GoogleTest)

add_executable(TableKernelGPUTest TableKernelGPUTest.cc)
target_link_libraries(TableKernelGPUTest gtest gtest_main)

add_test(
    NAME TableKernelGPUTest
    COMMAND TableKernelGPUTest
)
