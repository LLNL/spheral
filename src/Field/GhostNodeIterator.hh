//---------------------------------Spheral++----------------------------------//
// GhostNodeIterator -- The version of the NodeIterator that goes over all
// ghost nodes in a list of NodeLists.
//
// Created by J. Michael Owen, Mon Mar 17 22:56:00 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_GhostNodeIterator_hh__
#define __Spheral_GhostNodeIterator_hh__

#include "NodeIteratorBase.hh"

#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;

template<typename Dimension>
class GhostNodeIterator: public NodeIteratorBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors and destructors.
  GhostNodeIterator();
  GhostNodeIterator(typename std::vector<NodeList<Dimension>*>::const_iterator iter,
                    typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                    typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                    int nodeID = 0);
  GhostNodeIterator(const GhostNodeIterator& rhs);

  // Destructor.
  virtual ~GhostNodeIterator();

  // Increment the iterator.
  GhostNodeIterator& operator++();

  // Extended valid test for InternalNodeIterators.
  virtual bool valid() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  using NodeIteratorBase<Dimension>::mNodeID;
  using NodeIteratorBase<Dimension>::mFieldID;
  using NodeIteratorBase<Dimension>::mNodeListBegin;
  using NodeIteratorBase<Dimension>::mNodeListEnd;
  using NodeIteratorBase<Dimension>::mNodeListItr;
};

}

#include "GhostNodeIteratorInline.hh"

#endif
