//------------------------------------------------------------------------------
// Provide explicit versions of the template FieldList indexing test functions.
//------------------------------------------------------------------------------

#include "CXXTests/testFieldListIndexing.hh"

namespace Spheral {
namespace Testing {

// AllNodeIterators.
inline
std::string
testIndexScalarFieldListByAllNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  return testIndexByAllNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByAllNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  return testIndexByAllNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByAllNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  return testIndexByAllNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByAllNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  return testIndexByAllNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByAllNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  return testIndexByAllNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByAllNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  return testIndexByAllNodeIterators(dataBase, fieldList);
}

// InternalNodeIterators.
inline
std::string
testIndexScalarFieldListByInternalNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  return testIndexByInternalNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByInternalNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  return testIndexByInternalNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByInternalNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  return testIndexByInternalNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByInternalNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  return testIndexByInternalNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByInternalNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  return testIndexByInternalNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByInternalNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  return testIndexByInternalNodeIterators(dataBase, fieldList);
}

// GhostNodeIterators.
inline
std::string
testIndexScalarFieldListByGhostNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  return testIndexByGhostNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByGhostNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  return testIndexByGhostNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByGhostNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  return testIndexByGhostNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByGhostNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  return testIndexByGhostNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByGhostNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  return testIndexByGhostNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByGhostNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  return testIndexByGhostNodeIterators(dataBase, fieldList);
}

// MasterNodeIterators.
inline
std::string
testIndexScalarFieldListByMasterNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  return testIndexByMasterNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByMasterNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  return testIndexByMasterNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByMasterNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  return testIndexByMasterNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByMasterNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  return testIndexByMasterNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByMasterNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  return testIndexByMasterNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByMasterNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  return testIndexByMasterNodeIterators(dataBase, fieldList);
}

// CoarseNodeIterators.
inline
std::string
testIndexScalarFieldListByCoarseNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  return testIndexByCoarseNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByCoarseNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  return testIndexByCoarseNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByCoarseNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  return testIndexByCoarseNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByCoarseNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  return testIndexByCoarseNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByCoarseNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  return testIndexByCoarseNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByCoarseNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  return testIndexByCoarseNodeIterators(dataBase, fieldList);
}

// RefineNodeIterators.
inline
std::string
testIndexScalarFieldListByRefineNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
  return testIndexByRefineNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByRefineNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
                                             const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
  return testIndexByRefineNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByRefineNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
  return testIndexByRefineNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByRefineNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
                                             const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
  return testIndexByRefineNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexScalarFieldListByRefineNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
  return testIndexByRefineNodeIterators(dataBase, fieldList);
}

inline
std::string
testIndexVectorFieldListByRefineNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
                                             const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
  return testIndexByRefineNodeIterators(dataBase, fieldList);
}

// // CoarseNodeIterators (cache).
// inline
// std::string
// testCacheIndexScalarFieldListByCoarseNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
//                                                      const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
//   return testCacheIndexByCoarseNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexVectorFieldListByCoarseNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
//                                                      const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
//   return testCacheIndexByCoarseNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexScalarFieldListByCoarseNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
//                                                      const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
//   return testCacheIndexByCoarseNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexVectorFieldListByCoarseNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
//                                                      const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
//   return testCacheIndexByCoarseNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexScalarFieldListByCoarseNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
//                                                      const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
//   return testCacheIndexByCoarseNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexVectorFieldListByCoarseNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
//                                                      const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
//   return testCacheIndexByCoarseNodeIterators(dataBase, fieldList);
// }

// // RefineNodeIterators (cache).
// inline
// std::string
// testCacheIndexScalarFieldListByRefineNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
//                                                      const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Scalar>& fieldList) {
//   return testCacheIndexByRefineNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexVectorFieldListByRefineNodeIterators1d(const DataBase< Spheral::Dim<1> >& dataBase,
//                                                      const FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>& fieldList) {
//   return testCacheIndexByRefineNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexScalarFieldListByRefineNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
//                                                      const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Scalar>& fieldList) {
//   return testCacheIndexByRefineNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexVectorFieldListByRefineNodeIterators2d(const DataBase< Spheral::Dim<2> >& dataBase,
//                                                      const FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>& fieldList) {
//   return testCacheIndexByRefineNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexScalarFieldListByRefineNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
//                                                      const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Scalar>& fieldList) {
//   return testCacheIndexByRefineNodeIterators(dataBase, fieldList);
// }

// inline
// std::string
// testCacheIndexVectorFieldListByRefineNodeIterators3d(const DataBase< Spheral::Dim<3> >& dataBase,
//                                                      const FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>& fieldList) {
//   return testCacheIndexByRefineNodeIterators(dataBase, fieldList);
// }

}
}
