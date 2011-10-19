//---------------------------------Spheral++----------------------------------//
// NodeList -- An abstract base class for the NodeLists.
//
// We will define here the basic functionality we expect all NodeLists to 
// provide.
//
// Created by JMO, Wed Sep  8 21:54:50 PDT 1999
//----------------------------------------------------------------------------//

#include <algorithm>

#include "Geometry/Dimension.hh"
#include "NodeList.hh"
#include "NodeListRegistrar.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "Neighbor/Neighbor.hh"
#include "DataBase/State.hh"
#include "FileIO/FileIO.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/DBC.hh"
#include "Utilities/packElement.hh"

#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"

#include "TAU.h"

namespace Spheral {
namespace NodeSpace {

using namespace std;
using NeighborSpace::Neighbor;
using FileIOSpace::FileIO;
using KernelSpace::TableKernel;
using FieldSpace::FieldBase;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor with optional numInternal nodes, numGhost nodes, and name.
//------------------------------------------------------------------------------
template<typename Dimension>
NodeList<Dimension>::NodeList(std::string name,
                              const int numInternal,
                              const int numGhost,
                              const Scalar hmin,
                              const Scalar hmax,
                              const Scalar hminratio,
                              const Scalar nPerh,
                              const int maxNumNeighbors):
  mNumNodes(numInternal + numGhost),
  mFirstGhostNode(numInternal),
  mName(name),
  mMass(HydroFieldNames::mass),
  mPositions(HydroFieldNames::position),
  mVelocity(HydroFieldNames::velocity),
  mH(HydroFieldNames::H),
  mWork(HydroFieldNames::work),
  mhmin(hmin),
  mhmax(hmax),
  mhminratio(hminratio),
  mNodesPerSmoothingScale(nPerh),
  mMaxNumNeighbors(maxNumNeighbors),
  mFieldBaseList(),
  mNeighborPtr(0),
  mDummyList(),
  mRestart(DataOutput::registerWithRestart(*this, 10)) {
  REQUIRE(numInternal >= 0 && numGhost >= 0);
  mMass.setNodeList(*this);
  mPositions.setNodeList(*this);
  mVelocity.setNodeList(*this);
  mH.setNodeList(*this);
  mWork.setNodeList(*this);
  mDummyList.push_back(this);
  // It's never valid to have zero H's.
  mH = SymTensor::one;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NodeList<Dimension>::~NodeList() {
  // Loop over all the fields defined on this mesh, and destroy them.
  vector<FieldBase<Dimension>*> fieldBaseListCopy(mFieldBaseList);
  for (FieldBaseIterator fieldItr = fieldBaseListCopy.begin();
       fieldItr < fieldBaseListCopy.end(); ++fieldItr) {
    (*fieldItr)->unregisterNodeList();
  }

  // After we're done, all the field should have unregistered themselves
  // from the Node List.
  ENSURE(numFields() == 0);

}

//------------------------------------------------------------------------------
// Set the number of internal nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::numInternalNodes(int size) {
  CHECK(size >= 0);
  int numGhost = numGhostNodes();
  int oldFirstGhostNode = mFirstGhostNode;
  mFirstGhostNode = size;
  mNumNodes = size + numGhost;
  for (FieldBaseIterator fieldPtrItr = mFieldBaseList.begin();
       fieldPtrItr < mFieldBaseList.end(); ++fieldPtrItr) {
    (*fieldPtrItr)->resizeFieldInternal(size, oldFirstGhostNode);
  }
}

//------------------------------------------------------------------------------
// Set the number of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::numGhostNodes(int size) {
  CHECK(size >= 0);
  int numInternal = numInternalNodes();
  mNumNodes = numInternal + size;
  for (FieldBaseIterator fieldPtrItr = mFieldBaseList.begin();
       fieldPtrItr < mFieldBaseList.end(); ++fieldPtrItr) {
    (*fieldPtrItr)->resizeFieldGhost(size);
  }
}

//------------------------------------------------------------------------------
// Provide NodeIterators for all nodes in this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
AllNodeIterator<Dimension>
NodeList<Dimension>::nodeBegin() const {
  return AllNodeIterator<Dimension>(mDummyList.begin(),
                                    mDummyList.begin(),
                                    mDummyList.end());
}

template<typename Dimension>
AllNodeIterator<Dimension>
NodeList<Dimension>::nodeEnd() const {
  return AllNodeIterator<Dimension>(mDummyList.end(),
                                    mDummyList.begin(),
                                    mDummyList.end());
}

//------------------------------------------------------------------------------
// Node iterators for internal nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
InternalNodeIterator<Dimension>
NodeList<Dimension>::internalNodeBegin() const {
  return InternalNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end());
}

template<typename Dimension>
InternalNodeIterator<Dimension>
NodeList<Dimension>::internalNodeEnd() const {
  return InternalNodeIterator<Dimension>(mDummyList.end(),
                                         mDummyList.begin(),
                                         mDummyList.end());
}

//------------------------------------------------------------------------------
// Node iterators for ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
GhostNodeIterator<Dimension>
NodeList<Dimension>::ghostNodeBegin() const {
  return GhostNodeIterator<Dimension>(mDummyList.begin(),
                                      mDummyList.begin(),
                                      mDummyList.end(),
                                      firstGhostNode());
}

template<typename Dimension>
GhostNodeIterator<Dimension>
NodeList<Dimension>::ghostNodeEnd() const {
  return GhostNodeIterator<Dimension>(mDummyList.end(),
                                      mDummyList.begin(),
                                      mDummyList.end());
}

//------------------------------------------------------------------------------
// Node iterators for master neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
MasterNodeIterator<Dimension>
NodeList<Dimension>::masterNodeBegin() const {
  CHECK(mNeighborPtr != 0);
  if (mNeighborPtr->numMaster() > 0) {
    return MasterNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         mNeighborPtr->masterBegin());
  } else {
    return MasterNodeIterator<Dimension>(mDummyList.end(),
                                         mDummyList.begin(),
                                         mDummyList.end());
  }    
}

template<typename Dimension>
MasterNodeIterator<Dimension>
NodeList<Dimension>::masterNodeEnd() const {
  return MasterNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end());
}

//------------------------------------------------------------------------------
// Node iterators for coarse neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
CoarseNodeIterator<Dimension>
NodeList<Dimension>::coarseNodeBegin() const {
  CHECK(mNeighborPtr != 0);
  if (mNeighborPtr->numCoarse() > 0) {
    return CoarseNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         mNeighborPtr->coarseNeighborBegin());
  } else {
    return CoarseNodeIterator<Dimension>(mDummyList.end(),
                                         mDummyList.begin(),
                                         mDummyList.end());
  }    
}

template<typename Dimension>
CoarseNodeIterator<Dimension>
NodeList<Dimension>::coarseNodeEnd() const {
  return CoarseNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end());
}

//------------------------------------------------------------------------------
// Node iterators for refine neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
RefineNodeIterator<Dimension>
NodeList<Dimension>::refineNodeBegin() const {
  CHECK(mNeighborPtr != 0);
  if (mNeighborPtr->numRefine() > 0) {
    return RefineNodeIterator<Dimension>(mDummyList.begin(),
                                         mDummyList.begin(),
                                         mDummyList.end(),
                                         mNeighborPtr->refineNeighborBegin());
  } else {
    return RefineNodeIterator<Dimension>(mDummyList.end(),
                                         mDummyList.begin(),
                                         mDummyList.end());
  }    
}

template<typename Dimension>
RefineNodeIterator<Dimension>
NodeList<Dimension>::refineNodeEnd() const {
  return RefineNodeIterator<Dimension>(mDummyList.end(),
                                       mDummyList.begin(),
                                       mDummyList.end());
}


//------------------------------------------------------------------------------
// Set the mass field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
mass(const Field<Dimension, typename Dimension::Scalar>& m) {
  mMass = m;
  mMass.name(HydroFieldNames::mass);
}

//------------------------------------------------------------------------------
// Set the position field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
positions(const Field<Dimension, typename Dimension::Vector>& r) {
  mPositions = r;
  mPositions.name(HydroFieldNames::position);
}

//------------------------------------------------------------------------------
// Set the velocity field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
velocity(const Field<Dimension, typename Dimension::Vector>& v) {
  mVelocity = v;
  mVelocity.name(HydroFieldNames::velocity);
}

//------------------------------------------------------------------------------
// Set the H field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
Hfield(const Field<Dimension, typename Dimension::SymTensor>& H) {
  mH = H;
  mH.name(HydroFieldNames::H);
}

//------------------------------------------------------------------------------
// Set the work field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
work(const Field<Dimension, typename Dimension::Scalar>& w) {
  mWork = w;
  mWork.name(HydroFieldNames::work);
}

//------------------------------------------------------------------------------
// Compute and return the inverse field of H tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
Hinverse(Field<Dimension, typename Dimension::SymTensor>& field) const {
  REQUIRE(field.nodeListPtr() == this);
  for (int i = 0; i < numInternalNodes(); ++i) field(i) = mH(i).Inverse();
  field.name("H inverse");
}

//------------------------------------------------------------------------------
// Return the number of fields registered on this nodelist.
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeList<Dimension>::numFields() const {
  return mFieldBaseList.size();
}

//------------------------------------------------------------------------------
// Register a field with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::registerField(FieldBase<Dimension>& field) const {
  if (haveField(field)) {
    cerr << "WARNING: Attempt to register field " << &field
         << " with NodeList " << this << " that already has it." 
         << endl;
  } else {
    CHECK(&field);
    CHECK(&mFieldBaseList);
    mFieldBaseList.push_back(&field);
  }
}

//------------------------------------------------------------------------------
// Unregister a field that is listed with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::unregisterField(FieldBase<Dimension>& field) const {
  if (!haveField(field)) {
    cerr << "WARNING: Attempt to unregister field " << &field
         << " from a NodeList " << this << " that does not recognize it." 
         << endl;
  } else {
    FieldBaseIterator fieldPtrItr = find(mFieldBaseList.begin(),
                                         mFieldBaseList.end(),
                                         &field);
    mFieldBaseList.erase(fieldPtrItr);
  }
}

//------------------------------------------------------------------------------
// Check if the given field is registered with this NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NodeList<Dimension>::haveField(const FieldBase<Dimension>& field) const {
  const_FieldBaseIterator fieldPtrItr = find(mFieldBaseList.begin(),
                                             mFieldBaseList.end(),
                                             &field);
  return fieldPtrItr != mFieldBaseList.end();
}

//------------------------------------------------------------------------------
// Return the type of node (internal or ghost).
//------------------------------------------------------------------------------
template<typename Dimension>
NodeType
NodeList<Dimension>::nodeType(int i) const {
  CHECK(i >=0 && i < numNodes());
  if (i < firstGhostNode()) {
    return InternalNode;
  } else {
    return GhostNode;
  }
}

//------------------------------------------------------------------------------
// The index of the first ghost node.
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeList<Dimension>::firstGhostNode() const {
  CHECK(mFirstGhostNode >=0 && mFirstGhostNode <= numNodes());
  return mFirstGhostNode;
}

//------------------------------------------------------------------------------
// Access the neighbor object.
//------------------------------------------------------------------------------
template<typename Dimension>
Neighbor<Dimension>&
NodeList<Dimension>::neighbor() const {
  CHECK(mNeighborPtr != 0);
  return *mNeighborPtr;
}

//------------------------------------------------------------------------------
// Register the given neighbor object with this node list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
registerNeighbor(Neighbor<Dimension>& neighbor) {
  CHECK(&neighbor);
  mNeighborPtr = &neighbor;
}

//------------------------------------------------------------------------------
// Unregister the current neighbor object from this node list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::unregisterNeighbor() {
  mNeighborPtr = 0;
}

//------------------------------------------------------------------------------
// Delete the given node indicies from the NodeList (including all Fields).
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
deleteNodes(const vector<int>& nodeIDs) {

  // TAU timers.
  TAU_PROFILE("NodeList::", "deleteNodes(IDs)", TAU_USER);

  // First sort and make sure all node IDs are valid.
  vector<int> uniqueIDs(nodeIDs);
  sort(uniqueIDs.begin(), uniqueIDs.end());
  vector<int>::iterator uniqueEnd = unique(uniqueIDs.begin(), uniqueIDs.end());
  uniqueIDs.erase(uniqueEnd, uniqueIDs.end());
  CHECK(uniqueIDs.size() <= numNodes());
  if (uniqueIDs.size() > 0) 
    CHECK(uniqueIDs[0] >= 0 && uniqueIDs.back() < this->numNodes());

  // Determine how many internal, ghost, and total nodes we should end with.
  vector<int>::iterator ghostDeleteItr = uniqueIDs.begin();
  while (ghostDeleteItr < uniqueIDs.end() &&
         *ghostDeleteItr < mFirstGhostNode) ++ghostDeleteItr;
  CHECK(ghostDeleteItr >= uniqueIDs.begin() && ghostDeleteItr <= uniqueIDs.end());
  const int numInternalNodesRemoved = distance(uniqueIDs.begin(), ghostDeleteItr);
  CHECK(numInternalNodesRemoved <= numInternalNodes());
  mNumNodes -= uniqueIDs.size();
  mFirstGhostNode -= numInternalNodesRemoved;
  CHECK(mNumNodes >= 0);
  CHECK(mFirstGhostNode >= 0 && mFirstGhostNode <= mNumNodes);

  // Now iterate over the Fields defined on this NodeList, and remove the appropriate
  // elements from each.
  for (typename vector<FieldBase<Dimension>*>::iterator fieldItr = mFieldBaseList.begin();
       fieldItr != mFieldBaseList.end();
       ++fieldItr) (*fieldItr)->deleteElements(uniqueIDs);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  for (typename vector<FieldBase<Dimension>*>::iterator fieldItr = mFieldBaseList.begin();
       fieldItr < mFieldBaseList.end();
       ++fieldItr) {
    ENSURE((*fieldItr)->size() == mNumNodes);
  }
  END_CONTRACT_SCOPE;
  
}

//------------------------------------------------------------------------------
// Pack the Field values for the given node indicies as appends onto the given
// vector of Scalars.
//------------------------------------------------------------------------------
template<typename Dimension>
list< vector<char> >
NodeList<Dimension>::
packNodeFieldValues(const vector<int>& nodeIDs) const {

  // TAU timers.
  TAU_PROFILE("NodeList::", "packNodeFieldValues(IDs)", TAU_USER);

  // Prepare the result.
  list< vector<char> > result;

  // Sort and make sure all node IDs are valid.
  vector<int> uniqueIDs(nodeIDs);
  sort(uniqueIDs.begin(), uniqueIDs.end());
  vector<int>::iterator uniqueEnd = unique(uniqueIDs.begin(), uniqueIDs.end());
  uniqueIDs.erase(uniqueEnd, uniqueIDs.end());
  CHECK(uniqueIDs.size() <= numNodes());
  if (uniqueIDs.size() > 0) 
    CHECK(uniqueIDs[0] >= 0 && uniqueIDs.back() <= this->numNodes());

  // Iterate over all the Fields defined on this NodeList, and append it's packed 
  // field values to the stack.
  for (typename vector<FieldBase<Dimension>*>::const_iterator fieldItr = mFieldBaseList.begin();
       fieldItr != mFieldBaseList.end();
       ++fieldItr) {
    result.push_back((*fieldItr)->packValues(uniqueIDs));
  }

  ENSURE(result.size() == mFieldBaseList.size());
  return result;
}

//------------------------------------------------------------------------------
// Append the requested number of nodes onto this FieldList, and fill in the
// Field values using the provided packed data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
appendInternalNodes(const int numNewNodes,
                    const list< vector<char> >& packedFieldValues) {

  // TAU timers.
  TAU_PROFILE("NodeList::", "appendInternalValues(n, bufs)", TAU_USER);

  REQUIRE(packedFieldValues.size() == mFieldBaseList.size());

  // We only work if there are new nodes.
  if (numNewNodes > 0) {

    // Begin by resizing this NodeList appropriately.
    const int beginInsertionIndex = numInternalNodes();
    numInternalNodes(beginInsertionIndex + numNewNodes);
    CHECK(numInternalNodes() == beginInsertionIndex + numNewNodes);

    // Loop over each Field, and have them fill in the new values from the
    // packed char buffers.
    typename list< vector<char> >::const_iterator bufItr = packedFieldValues.begin();
    for (typename vector<FieldBase<Dimension>*>::iterator fieldItr = mFieldBaseList.begin();
         fieldItr != mFieldBaseList.end();
         ++fieldItr, ++bufItr) {
      CHECK(bufItr != packedFieldValues.end());
      (*fieldItr)->unpackValues(numNewNodes, 
                                beginInsertionIndex,
                                *bufItr);
    }

    // That's it.
    ENSURE(bufItr == packedFieldValues.end());
  }
}

//------------------------------------------------------------------------------
// Reorder the nodes according to the requested ordering.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
reorderNodes(const vector<int>& newOrdering) {

  // TAU timers.
  TAU_PROFILE("NodeList::", "reorderNodes(order)", TAU_USER);

  // The number of internal nodes.
  const int n = this->numInternalNodes();

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(newOrdering.size() == n);
    vector<int> tmp(newOrdering);
    sort(tmp.begin(), tmp.end());
    for (int i = 0; i != n; ++i) REQUIRE(tmp[i] == i);
  }
  END_CONTRACT_SCOPE;

  // Make sure we're not carting around ghost nodes.
  this->numGhostNodes(0);

  // Pack up all the current nodal field values.
  list<vector<char> > packedFieldValues;
  for (typename vector<FieldBase<Dimension>*>::const_iterator fieldItr = mFieldBaseList.begin();
       fieldItr != mFieldBaseList.end();
       ++fieldItr) packedFieldValues.push_back((*fieldItr)->packValues(newOrdering));
  CHECK(packedFieldValues.size() == mFieldBaseList.size());

  // Zap out all the current nodes.
  this->numInternalNodes(0);

  // Now unpack in the desired order.
  this->appendInternalNodes(n, packedFieldValues);

  // Post-conditions.
  ENSURE(this->numInternalNodes() == n);
}

//------------------------------------------------------------------------------
// Dump the current state of the NodeList to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Dump the name of the NodeList.
  file.write(name(), pathName + "/name");

  // Dump the number of internal nodes (we assume that ghost information
  // does not need to be stored).
  file.write(numInternalNodes(), pathName + "/numNodes");

  // Dump each of the internal fields of the NodeList.
  file.write(mMass, pathName + "/mass");
  file.write(mPositions, pathName + "/positions");
  file.write(mVelocity, pathName + "/velocity");
  file.write(mH, pathName + "/H");
  file.write(mWork, pathName + "/work");
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeList<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // Restore the name of the NodeList.
  file.read(mName, pathName + "/name");

  // Read and reset the number of internal nodes.
  int numNodes;
  file.read(numNodes, pathName + "/numNodes");
  numInternalNodes(numNodes);

  // Now we can restore each of the internal fields of the NodeList.
  file.read(mMass, pathName + "/mass");
  file.read(mPositions, pathName + "/positions");
  file.read(mVelocity, pathName + "/velocity");
  file.read(mH, pathName + "/H");
  file.read(mWork, pathName + "/work");

  // The neighbor object doesn't actually write out state, but does need to be
  // reinitialized with the new NodeList state.
  mNeighborPtr->updateNodes();
}

// //------------------------------------------------------------------------------
// // Notify all Fields registered on this NodeList to cache their coarse neighbor
// // values.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// NodeList<Dimension>::
// notifyFieldsCacheCoarseValues() const {
//   for (typename vector<FieldBase<Dimension>*>::iterator fieldItr = mFieldBaseList.begin();
//        fieldItr < mFieldBaseList.end();
//        ++fieldItr) {
//     (*fieldItr)->notifyNewCoarseNodes();
//   }
// }

// //------------------------------------------------------------------------------
// // Notify all Fields registered on this NodeList to cache their refine neighbor
// // values.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// NodeList<Dimension>::
// notifyFieldsCacheRefineValues() const {
//   for (typename vector<FieldBase<Dimension>*>::iterator fieldItr = mFieldBaseList.begin();
//        fieldItr < mFieldBaseList.end();
//        ++fieldItr) {
//     (*fieldItr)->notifyNewRefineNodes();
//   }
// }

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class NodeList< Dim<1> >;
template class NodeList< Dim<2> >;
template class NodeList< Dim<3> >;
}
}

