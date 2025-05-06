spheral_add_test(
  NAME managed_vector_perf
  SOURCES managed_vector_perf.cc
  #DEPENDS_ON Spheral_Geometry Spheral_Utilities Spheral_Distributed
  DEBUG_LINKER
)

spheral_add_test(
  NAME managed_vector_tests
  SOURCES managed_vector_tests.cc
  DEPENDS_ON Spheral_Geometry Spheral_Utilities Spheral_Distributed
  DEBUG_LINKER
)
