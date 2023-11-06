
add_subdirectory(DeviceTestLib)

if(ENABLE_TESTS)
  if(ENABLE_DEV_BUILD)
    get_property(SPHERAL_DEPENDS GLOBAL PROPERTY SPHERAL_OBJ_LIBS)
  else()
    set(SPHERAL_DEPENDS Spheral_CXX)
  endif()
  blt_add_executable( NAME spheral_cuda_test
                      SOURCES Spheral_CUDA_Test.cc
                      DEPENDS_ON ${SPHERAL_DEPENDS} ${spheral_blt_cxx_depends})
endif()
