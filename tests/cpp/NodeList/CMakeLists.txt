spheral_add_test(
  NAME nodelistview_tests
  SOURCES nodelistview_tests.cc
  DEPENDS_ON Spheral_Distributed Spheral_NodeList Spheral_Field Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_Utilities Spheral_FileIO
  #DEBUG_LINKER
)
