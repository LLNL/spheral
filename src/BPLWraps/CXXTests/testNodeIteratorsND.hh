//------------------------------------------------------------------------------
// Provide explicit versions of the template node iterator test functions.
//------------------------------------------------------------------------------

#include <string>
#include "CXXTests/testNodeIterators.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace Testing {

using DataBaseSpace::DataBase;

// AllNodeIterators.
inline
std::string
testGlobalAllNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase) {
  return testGlobalAllNodeIterators(dataBase);
}

inline
std::string
testGlobalAllNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase) {
  return testGlobalAllNodeIterators(dataBase);
}

inline
std::string
testGlobalAllNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase) {
  return testGlobalAllNodeIterators(dataBase);
}

// InternalNodeIterators.
inline
std::string
testGlobalInternalNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase) {
  return testGlobalInternalNodeIterators(dataBase);
}

inline
std::string
testGlobalInternalNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase) {
  return testGlobalInternalNodeIterators(dataBase);
}

inline
std::string
testGlobalInternalNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase) {
  return testGlobalInternalNodeIterators(dataBase);
}

// GhostNodeIterators.
inline
std::string
testGlobalGhostNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase) {
  return testGlobalGhostNodeIterators(dataBase);
}

inline
std::string
testGlobalGhostNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase) {
  return testGlobalGhostNodeIterators(dataBase);
}

inline
std::string
testGlobalGhostNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase) {
  return testGlobalGhostNodeIterators(dataBase);
}

// MasterNodeIterators.
inline
std::string
testGlobalMasterNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase) {
  return testGlobalMasterNodeIterators(dataBase);
}

inline
std::string
testGlobalMasterNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase) {
  return testGlobalMasterNodeIterators(dataBase);
}

inline
std::string
testGlobalMasterNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase) {
  return testGlobalMasterNodeIterators(dataBase);
}

// CoarseNodeIterators.
inline
std::string
testGlobalCoarseNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase) {
  return testGlobalCoarseNodeIterators(dataBase);
}

inline
std::string
testGlobalCoarseNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase) {
  return testGlobalCoarseNodeIterators(dataBase);
}

inline
std::string
testGlobalCoarseNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase) {
  return testGlobalCoarseNodeIterators(dataBase);
}

// RefineNodeIterators.
inline
std::string
testGlobalRefineNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase) {
  return testGlobalRefineNodeIterators(dataBase);
}

inline
std::string
testGlobalRefineNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase) {
  return testGlobalRefineNodeIterators(dataBase);
}

inline
std::string
testGlobalRefineNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase) {
  return testGlobalRefineNodeIterators(dataBase);
}

}
}
