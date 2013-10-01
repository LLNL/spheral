//------------------------------------------------------------------------------
// testNodeIterators
// A collection of test functions for the node iterators.
//
// Created by JMO, Thu Mar 25 14:21:41 2004
//------------------------------------------------------------------------------
#ifndef __Spheral_testNodeIterators_hh__
#define __Spheral_testNodeIterators_hh__

#ifndef __GCCXML__
#include <string>
#else
#include "Utilities/fakestl.hh"
#endif

namespace Spheral {
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
}

namespace Spheral {
namespace Testing {

//------------------------------------------------------------------------------
// Test global AllNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalAllNodeIterators(const DataBaseSpace::DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global InternalNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalInternalNodeIterators(const DataBaseSpace::DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global GhostNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalGhostNodeIterators(const DataBaseSpace::DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global MasterNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalMasterNodeIterators(const DataBaseSpace::DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global CoarseNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalCoarseNodeIterators(const DataBaseSpace::DataBase<Dimension>& dataBase);

//------------------------------------------------------------------------------
// Test global RefineNodeIterators.
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
testGlobalRefineNodeIterators(const DataBaseSpace::DataBase<Dimension>& dataBase);

}
}

#endif
