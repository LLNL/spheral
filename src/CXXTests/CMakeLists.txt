get_property(SPHERAL_CXX_LIBS GLOBAL PROPERTY SPHERAL_CXX_LIBS)

add_subdirectory(DeviceTestLib)

blt_add_executable( NAME spheral_cuda_test
                    SOURCES Spheral_CUDA_Test.cc
                    DEPENDS_ON Spheral_CXX ${SPHERAL_CXX_DEPENDS})
