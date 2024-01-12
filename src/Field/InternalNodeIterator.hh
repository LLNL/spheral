//---------------------------------Spheral++----------------------------------//
// InternalNodeIterator -- The version of the NodeIterator that goes over all
// internal nodes in a list of NodeLists.
//
// Created by J. Michael Owen, Mon Mar 17 23:13:08 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_InternalNodeIterator_hh__
#define __Spheral_InternalNodeIterator_hh__

#include "NodeIteratorBase.hh"

#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;

template<typename Dimension>
class InternalNodeIterator: public NodeIteratorBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors and destructors.
  InternalNodeIterator();
  InternalNodeIterator(typename std::vector<NodeList<Dimension>*>::const_iterator iter,
                       typename std::vector<NodeList<Dimension>*>::const_iterator nodeListBegin,
                       typename std::vector<NodeList<Dimension>*>::const_iterator nodeListEnd,
                       int nodeID = 0);
  InternalNodeIterator(const InternalNodeIterator& rhs);

  // Destructor.
  virtual ~InternalNodeIterator();

  // Increment the iterator.
  InternalNodeIterator& operator++();

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

#include "InternalNodeIteratorInline.hh"

#endif
