
add_subdirectory(DeviceTestLib)

if(ENABLE_TESTS)
  blt_add_executable( NAME spheral_cuda_test
                      SOURCES Spheral_CUDA_Test.cc
                      DEPENDS_ON Spheral_DeviceTestLib ${spheral_blt_cxx_depends})
endif()
