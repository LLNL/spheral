spheral_add_test(
  NAME managed_vector_tests
  SOURCES managed_vector_tests.cc
)

spheral_add_test(
  NAME field_parallel_inheritance_tests
  SOURCES field_parallel_inheritance_tests.cc
)

#spheral_add_test(
#  NAME field_tests
#  SOURCES field_tests.cc
#  DEPENDS_ON Spheral_NodeList Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_Utilities
#  #DEBUG_LINKER
#)
#
#spheral_add_test(
#  NAME fieldview_tests
#  SOURCES fieldview_tests.cc
#  #DEBUG_LINKER
#)
