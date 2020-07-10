//------------------------------------------------------------------------------
// testFieldListIndexing
// A collection of test functions for indexing into FieldLists with 
// NodeIterators.
//
// Created by JMO, Wed Apr 14 11:10:22 2004
//------------------------------------------------------------------------------
#ifndef __Spheral_testFieldListIndexing_hh__
#define __Spheral_testFieldListIndexing_hh__

#include <map>
#include <string>
#include <sstream>

#include "Field/NodeIterators.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Workhorse method to test that the indexed FieldList values for a given range
// of NodeIterators matches the expected field values.
//------------------------------------------------------------------------------
template<typename FieldListType, typename IteratorType>
std::string
testFieldListIndexing(IteratorType beginIterator,
                      IteratorType endIterator,
                      const FieldListType& fieldList) {

  // Iterate over the range of the given iterator.
  for (IteratorType itr = beginIterator; itr != endIterator; ++itr) {

    // Check that index operator for the iterator returns the expected
    // sub field value.
    if (fieldList(itr) != (*fieldList[itr.fieldID()])[itr.nodeID()]) {
      std::stringstream message;
      message << "Wrong field value for "
              << itr.fieldID() << " "
              << itr.nodeID() << " : "
              << fieldList(itr) << " != "
              << (*fieldList[itr.fieldID()])[itr.nodeID()] << std::endl;
      return message.str();
    }

  }

  return "OK";
}

// //------------------------------------------------------------------------------
// // Similar method testing the cached coarse FieldList values for a given range
// // of CoarseNodeIterators.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// std::string
// testCoarseCacheFieldListIndexing(CoarseNodeIterator<Dimension> beginIterator,
//                                  CoarseNodeIterator<Dimension> endIterator,
//                                  const FieldList<Dimension, DataType>& fieldList) {

//   // Make sure the field list has cached the coarse values.
//   fieldList.cacheCoarseValues();

//   // Iterate over the range of the given iterator.
//   typedef CoarseNodeIterator<Dimension> IteratorType;
//   for (IteratorType itr = beginIterator; itr != endIterator; ++itr) {

//     // Check that index operator for the iterator returns the expected
//     // sub field value.
//     if (fieldList.coarseCache(itr) != (*fieldList[itr.fieldID()])[itr.nodeID()]) {
//       std::stringstream message;
//       message << "Wrong field value for "
//               << itr.fieldID() << " "
//               << itr.nodeID() << " : "
//               << fieldList.coarseCache(itr) << " != "
//               << (*fieldList[itr.fieldID()])[itr.nodeID()] << std::endl;
//       return message.str();
//     }

//   }

//   return "OK";
// }

// //------------------------------------------------------------------------------
// // Similar method testing the cached refine FieldList values for a given range
// // of RefineNodeIterators.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// std::string
// testRefineCacheFieldListIndexing(RefineNodeIterator<Dimension> beginIterator,
//                                  RefineNodeIterator<Dimension> endIterator,
//                                  const FieldList<Dimension, DataType>& fieldList) {

//   // Make sure the field list has cached the refine values.
//   fieldList.cacheRefineValues();

//   // Iterate over the range of the given iterator.
//   typedef RefineNodeIterator<Dimension> IteratorType;
//   for (IteratorType itr = beginIterator; itr != endIterator; ++itr) {

//     // Check that index operator for the iterator returns the expected
//     // sub field value.
//     if (fieldList.refineCache(itr) != (*fieldList[itr.fieldID()])[itr.nodeID()]) {
//       std::stringstream message;
//       message << "Wrong field value for "
//               << itr.fieldID() << " "
//               << itr.nodeID() << " : "
//               << fieldList.refineCache(itr) << " != "
//               << (*fieldList[itr.fieldID()])[itr.nodeID()] << std::endl;
//       return message.str();
//     }

//   }

//   return "OK";
// }

//------------------------------------------------------------------------------
// Test indexing by global AllNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::string
testIndexByAllNodeIterators(const DataBase<Dimension>& dataBase,
                            const FieldList<Dimension, DataType>& fieldList) {
  return testFieldListIndexing(dataBase.nodeBegin(),
                               dataBase.nodeEnd(),
                               fieldList);
}

//------------------------------------------------------------------------------
// Test indexing by global InternalNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::string
testIndexByInternalNodeIterators(const DataBase<Dimension>& dataBase,
                                 const FieldList<Dimension, DataType>& fieldList) {
  return testFieldListIndexing(dataBase.internalNodeBegin(),
                               dataBase.internalNodeEnd(),
                               fieldList);
}

//------------------------------------------------------------------------------
// Test indexing by global GhostNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::string
testIndexByGhostNodeIterators(const DataBase<Dimension>& dataBase,
                              const FieldList<Dimension, DataType>& fieldList) {
  return testFieldListIndexing(dataBase.ghostNodeBegin(),
                               dataBase.ghostNodeEnd(),
                               fieldList);
}

//------------------------------------------------------------------------------
// Test indexing by global MasterNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::string
testIndexByMasterNodeIterators(const DataBase<Dimension>& dataBase,
                               const FieldList<Dimension, DataType>& fieldList) {
  return testFieldListIndexing(dataBase.masterNodeBegin(),
                               dataBase.masterNodeEnd(),
                               fieldList);
}

//------------------------------------------------------------------------------
// Test indexing by global CoarseNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::string
testIndexByCoarseNodeIterators(const DataBase<Dimension>& dataBase,
                              const FieldList<Dimension, DataType>& fieldList) {
  return testFieldListIndexing(dataBase.coarseNodeBegin(),
                               dataBase.coarseNodeEnd(),
                               fieldList);
}

// template<typename Dimension, typename DataType>
// inline
// std::string
// testCacheIndexByCoarseNodeIterators(const DataBase<Dimension>& dataBase,
//                                     const FieldList<Dimension, DataType>& fieldList) {
//   return testCoarseCacheFieldListIndexing(dataBase.coarseNodeBegin(),
//                                           dataBase.coarseNodeEnd(),
//                                           fieldList);
// }

//------------------------------------------------------------------------------
// Test indexing by global RefineNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::string
testIndexByRefineNodeIterators(const DataBase<Dimension>& dataBase,
                               const FieldList<Dimension, DataType>& fieldList) {
  return testFieldListIndexing(dataBase.refineNodeBegin(),
                               dataBase.refineNodeEnd(),
                               fieldList);
}

// template<typename Dimension, typename DataType>
// inline
// std::string
// testCacheIndexByRefineNodeIterators(const DataBase<Dimension>& dataBase,
//                                     const FieldList<Dimension, DataType>& fieldList) {
//   return testRefineCacheFieldListIndexing(dataBase.refineNodeBegin(),
//                                           dataBase.refineNodeEnd(),
//                                           fieldList);
// }

}

#endif
