//---------------------------------Spheral++----------------------------------//
// CoarseNodeIterator -- The version of the NodeIterator that goes over all
// coarse nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Tue Mar 18 21:32:13 PST 2003
//----------------------------------------------------------------------------//
#ifndef CoarseNodeIterator_HH
#define CoarseNodeIterator_HH

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
class CoarseNodeIterator: public NodeIteratorBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors and destructors.
  CoarseNodeIterator();
#ifndef __GCCXML__
  CoarseNodeIterator(typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator iter,
                     typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListEnd);
  CoarseNodeIterator(typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator iter,
                     typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListEnd,
                     std::vector<int>::const_iterator IDItr);
#endif
  CoarseNodeIterator(const CoarseNodeIterator& rhs);

  // Destructor.
  virtual ~CoarseNodeIterator();

  // Increment the iterator.
  CoarseNodeIterator& operator++();

  // Valid test.
  virtual bool valid() const;

private:
  //---------------------------- Private Interface ----------------------------//
  // Internal method to initialize the state.
  void initialize(typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListItr,
                  typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListBegin,
                  typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListEnd,
                  std::vector<int>::const_iterator IDItr);

#ifndef __GCCXML__
  // The iterator to the current coarse node ID.
  typename std::vector<int>::const_iterator mNodeIDItr;

protected:
  //--------------------------- Protected Interface ---------------------------//
  using NodeIteratorBase<Dimension>::mNodeID;
  using NodeIteratorBase<Dimension>::mFieldID;
  using NodeIteratorBase<Dimension>::mNodeListBegin;
  using NodeIteratorBase<Dimension>::mNodeListEnd;
  using NodeIteratorBase<Dimension>::mNodeListItr;
#endif
};

}

#ifndef __GCCXML__
#include "CoarseNodeIteratorInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class CoarseNodeIterator;
}

#endif
