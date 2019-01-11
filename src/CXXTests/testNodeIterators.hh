//------------------------------------------------------------------------------
// testNodeIterators
// A collection of test functions for the node iterators.
//
// Created by JMO, Thu Mar 25 14:21:41 2004
//------------------------------------------------------------------------------
#ifndef __Spheral_testNodeIterators_hh__
#define __Spheral_testNodeIterators_hh__

#include <string>

namespace Spheral {
  template<typename Dimension> class DataBase;
}

namespace Spheral {

//------------------------------------------------------------------------------
// Test global AllNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalAllNodeIterators(const DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global InternalNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalInternalNodeIterators(const DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global GhostNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalGhostNodeIterators(const DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global MasterNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalMasterNodeIterators(const DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global CoarseNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalCoarseNodeIterators(const DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global RefineNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalRefineNodeIterators(const DataBase<Dimension>& dataBase);

}

#endif
