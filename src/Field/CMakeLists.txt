blt_add_library(
  NAME Field_spheral
  HEADERS AllNodeIterator.hh
	  CoarseNodeIterator.hh
	  Field.hh
	  FieldBase.hh
	  FieldList.hh
	  FieldListBase.hh
	  FieldListSet.hh
	  GhostNodeIterator.hh
	  InternalNodeIterator.hh
	  MasterNodeIterator.hh
	  NodeIteratorBase.hh
	  NodeIterators.hh
	  RefineNodeIterator.hh
  SOURCES AllNodeIterator.cc
	  CoarseNodeIterator.cc
	  FieldListBase.cc
	  GhostNodeIterator.cc
	  InternalNodeIterator.cc
	  MasterNodeIterator.cc
	  NodeIteratorBase.cc
	  RefineNodeIterator.cc
  OBJECT FALSE
  DEPENDS_ON Utilities_spheral)
