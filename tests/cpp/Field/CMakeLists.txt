spheral_add_test(
  NAME field_tests
  SOURCES field_tests.cc
  DEPENDS_ON Spheral_NodeList Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_Utilities
  #DEBUG_LINKER
)
