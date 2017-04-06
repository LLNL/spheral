//---------------------------------Spheral++----------------------------------//
// AllNodeIterator -- The version of the NodeIterator that goes over all nodes
// in a list of NodeLists.
//
// Created by J. Michael Owen, Mon Mar 17 17:03:38 PST 2003
//----------------------------------------------------------------------------//
#ifndef AllNodeIterator_HH
#define AllNodeIterator_HH

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "NodeIteratorBase.hh"

namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
}

namespace Spheral {

template<typename Dimension>
class AllNodeIterator: public NodeIteratorBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors and destructors.
  AllNodeIterator();
  AllNodeIterator(typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator itr,
                  typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListBegin,
                  typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListEnd,
                  int nodeID = 0);
  AllNodeIterator(const AllNodeIterator& rhs);

  // Destructor.
  virtual ~AllNodeIterator();

  // Increment the iterator.
  AllNodeIterator& operator++();

protected:
  //--------------------------- Protected Interface ---------------------------//
#ifndef __GCCXML__
  using NodeIteratorBase<Dimension>::mNodeID;
  using NodeIteratorBase<Dimension>::mFieldID;
  using NodeIteratorBase<Dimension>::mNodeListBegin;
  using NodeIteratorBase<Dimension>::mNodeListEnd;
  using NodeIteratorBase<Dimension>::mNodeListItr;
#endif
};

}

#ifndef __GCCXML__
#include "AllNodeIteratorInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class AllNodeIterator;
}

#endif
