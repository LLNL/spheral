//---------------------------------Spheral++----------------------------------//
// CoarseNodeIterator -- The version of the NodeIterator that goes over all
// coarse nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Tue Mar 18 21:32:13 PST 2003
//----------------------------------------------------------------------------//
#ifndef CoarseNodeIterator_HH
#define CoarseNodeIterator_HH

#include "NodeIteratorBase.hh"
#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;

template<typename Dimension>
class CoarseNodeIterator: public NodeIteratorBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors and destructors.
  CoarseNodeIterator();
  CoarseNodeIterator(typename std::vector<NodeList<Dimension>*>::const_iterator iter,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                     const std::vector<std::vector<int>>& coarseNeighbors);
  CoarseNodeIterator(typename std::vector<NodeList<Dimension>*>::const_iterator iter,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                     std::vector<int>::const_iterator IDItr,
                     const std::vector<std::vector<int>>& coarseNeighbors);
  CoarseNodeIterator(const CoarseNodeIterator& rhs);

  // Destructor.
  virtual ~CoarseNodeIterator();

  // Increment the iterator.
  CoarseNodeIterator& operator++();

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
                  const std::vector<std::vector<int>>& coarseNeighbors);

  // The iterator to the current coarse node ID.
  typename std::vector<int>::const_iterator mNodeIDItr;
  std::vector<std::vector<int>> mCoarseNeighbors;
};

}

#include "CoarseNodeIteratorInline.hh"

#endif
