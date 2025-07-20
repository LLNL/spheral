spheral_add_test(
  NAME geomvector_test
  SOURCES geomvector_tests.cc
  DEPENDS_ON Spheral_Utilities
  #DEBUG_LINKER
)

spheral_add_test(
  NAME rankntensor_test
  SOURCES rankntensor_tests.cc
  DEPENDS_ON Spheral_Utilities
  #DEBUG_LINKER
)
