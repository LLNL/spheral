
add_subdirectory(DeviceTestLib)

if(ENABLE_TESTS AND NOT ENABLE_DEV_BUILD)
  blt_add_executable( NAME spheral_cuda_test
                      SOURCES Spheral_CUDA_Test.cc
                      DEPENDS_ON Spheral_CXX ${spheral_blt_cxx_depends})
endif()
