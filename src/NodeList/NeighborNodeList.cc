//---------------------------------Spheral++----------------------------------//
// NeighborNodeList -- An abstract base class for the NodeLists to represent
//                  neighbors.
//
// Created by JMO, Sat Sep 18 10:50:42 PDT 1999
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "SmoothingScaleBase.hh"
#include "Material/EquationOfState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/safeInv.hh"
#include "NeighborNodeList.hh"

using std::vector;
using std::list;
using std::string;
using std::cerr;
using std::endl;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with given EOS object, along with optional numInternal nodes,
// numGhost nodes, and name.
//------------------------------------------------------------------------------
template<typename Dimension>
NeighborNodeList<Dimension>::
NeighborNodeList(string name,
                 const int numInternal,
                 const int numGhost,
                 const Scalar hmin,
                 const Scalar hmax,
                 const Scalar hminratio,
                 const Scalar nPerh,
                 const int maxNumNeighbors):
   NodeList<Dimension>(name, numInternal, numGhost, hmin, hmax, hminratio, nPerh),
   mMaxNumNeighbors(maxNumNeighbors),
   mNeighborPtr(0) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
NeighborNodeList<Dimension>::
~NeighborNodeList() {
}

//------------------------------------------------------------------------------
// Node iterators for master neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>
NodeList<Dimension>::masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const {
  REQUIRE(mNeighborPtr != 0);
  REQUIRE(masterLists.size() == 1);
  if (masterLists[0].size() > 0) {
    return MasterNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         masterLists[0].begin(),
                                         masterLists);
  } else {
    return this->masterNodeEnd();
  }
}

template<typename Dimension>
MasterNodeIterator<Dimension>
NodeList<Dimension>::masterNodeEnd() const {
  return MasterNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Node iterators for coarse neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>
NodeList<Dimension>::coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const {
  REQUIRE(mNeighborPtr != 0);
  REQUIRE(coarseNeighbors.size() == 1);
  if (coarseNeighbors[0].size() > 0) {
    return CoarseNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         coarseNeighbors[0].begin(),
                                         coarseNeighbors);
  } else {
    return this->coarseNodeEnd();
  }
}

template<typename Dimension>
CoarseNodeIterator<Dimension>
NodeList<Dimension>::coarseNodeEnd() const {
  return CoarseNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Node iterators for refine neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
RefineNodeIterator<Dimension>
NodeList<Dimension>::refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const {
  REQUIRE(mNeighborPtr != 0);
  REQUIRE(refineNeighbors.size() == 1);
  if (refineNeighbors[0].size() > 0) {
    return RefineNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         refineNeighbors[0].begin(),
                                         refineNeighbors);
  } else {
    return this->refineNodeEnd();
  }
}

template<typename Dimension>
RefineNodeIterator<Dimension>
NodeList<Dimension>::refineNodeEnd() const {
  return RefineNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Dump the current state of the NodeList to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NeighborNodeList<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  // Dump the ancestor class.
  NodeList<Dimension>::dumpState(file, pathName);
}

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NeighborNodeList<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  // Restore the ancestor class.
  NodeList<Dimension>::restoreState(file, pathName);
}  

}
