//---------------------------------Spheral++----------------------------------//
// MasterNodeIterator -- The version of the NodeIterator that goes over all
// master nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Mon Mar 17 23:13:08 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_MasterNodeIterator_hh__
#define __Spheral_MasterNodeIterator_hh__

#include "NodeIteratorBase.hh"
#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;

template<typename Dimension>
class MasterNodeIterator: public NodeIteratorBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors and destructors.
  MasterNodeIterator();
  MasterNodeIterator(typename std::vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                     const std::vector<std::vector<int>>& masterLists);
  MasterNodeIterator(typename std::vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                     std::vector<int>::const_iterator IDItr,
                     const std::vector<std::vector<int>>& masterLists);
  MasterNodeIterator(const MasterNodeIterator& rhs);

  // Destructor.
  virtual ~MasterNodeIterator();

  // Increment the iterator.
  MasterNodeIterator& operator++();

  // Valid test.
  virtual bool valid() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  using NodeIteratorBase<Dimension>::mNodeID;
  using NodeIteratorBase<Dimension>::mFieldID;
  using NodeIteratorBase<Dimension>::mNodeListBegin;
  using NodeIteratorBase<Dimension>::mNodeListEnd;
  using NodeIteratorBase<Dimension>::mNodeListItr;

private:
  //--------------------------- Private Interface ---------------------------//
  // Internal method to initialize the state.
  void initialize(typename std::vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                  typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                  typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                  std::vector<int>::const_iterator IDItr,
                  const std::vector<std::vector<int>>& masterLists);

  // The iterator to the current master node ID.
  typename std::vector<int>::const_iterator mNodeIDItr;
  std::vector<std::vector<int>> mMasterLists;
};

}

#include "MasterNodeIteratorInline.hh"

#endif
