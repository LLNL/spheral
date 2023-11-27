//---------------------------------Spheral++----------------------------------//
// RefineNodeIterator -- The version of the NodeIterator that goes over all
// refine nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Tue Mar 18 21:32:13 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_RefineNodeIterator_hh__
#define __Spheral_RefineNodeIterator_hh__

#include "NodeIteratorBase.hh"
#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;

template<typename Dimension>
class RefineNodeIterator: public NodeIteratorBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors and destructors.
  RefineNodeIterator();
  RefineNodeIterator(typename std::vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                     const std::vector<std::vector<int>>& refineNeighbors);
  RefineNodeIterator(typename std::vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                     std::vector<int>::const_iterator IDItr,
                     const std::vector<std::vector<int>>& refineNeighbors);
  RefineNodeIterator(const RefineNodeIterator& rhs);

  // Destructor.
  virtual ~RefineNodeIterator();

  // Increment the iterator.
  RefineNodeIterator& operator++();

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
  //---------------------------- Private Interface ----------------------------//
  // Internal method to initialize the state.
  void initialize(typename std::vector<NodeList<Dimension>*>::const_iterator nodeListItr,
                  typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                  typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                  std::vector<int>::const_iterator IDItr,
                  const std::vector<std::vector<int>>& refineNeighbors);

  // The iterator to the current refine node ID.
  typename std::vector<int>::const_iterator mNodeIDItr;
  std::vector<std::vector<int>> mRefineNeighbors;
};

}

#include "RefineNodeIteratorInline.hh"

#endif
