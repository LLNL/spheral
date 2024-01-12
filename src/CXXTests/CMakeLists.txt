
add_subdirectory(DeviceTestLib)

#if(ENABLE_TESTS)
  spheral_add_test( 
    NAME spheral_cuda_test
    SOURCES Spheral_CUDA_Test.cc
    DEPENDS_ON Spheral_Utilities Spheral_NodeList Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_DataBase Spheral_Neighbor
  )
#endif()
