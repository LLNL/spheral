spheral_add_test(
  NAME field_tests
  SOURCES field_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Distributed Spheral_NodeList Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_Utilities Spheral_FileIO
  #DEBUG_LINKER
)

spheral_add_test(
  NAME fieldview_tests
  SOURCES fieldview_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Distributed Spheral_NodeList Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_Utilities Spheral_FileIO
  #DEBUG_LINKER
)

spheral_add_test(
  NAME fieldview_datatype_tests
  SOURCES fieldview_datatype_tests.cc
  DEPENDS_ON Spheral_Field Spheral_NodeList Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_Utilities Spheral_Distributed Spheral_FileIO
  #DEBUG_LINKER
)

spheral_add_test(
  NAME fieldlist_tests
  SOURCES fieldlist_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Distributed Spheral_NodeList Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_Utilities Spheral_FileIO
  #DEBUG_LINKER
)

spheral_add_test(
  NAME fieldlistview_tests
  SOURCES fieldlistview_tests.cc
  DEPENDS_ON Spheral_Field Spheral_Distributed Spheral_NodeList Spheral_Geometry Spheral_Hydro Spheral_DataOutput Spheral_Utilities Spheral_FileIO
  #DEBUG_LINKER
)
