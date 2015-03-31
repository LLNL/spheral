//---------------------------------Spheral++----------------------------------//
// MasterNodeIterator -- The version of the NodeIterator that goes over all
// master nodes in a list of NodeLists based on their Neighbor states.
//
// Created by J. Michael Owen, Mon Mar 17 23:13:08 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_MasterNodeIterator_hh__
#define __Spheral_MasterNodeIterator_hh__

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
class MasterNodeIterator: public NodeIteratorBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors and destructors.
  MasterNodeIterator();
#ifndef __GCCXML__
  MasterNodeIterator(typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListItr,
                     typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListEnd);
  MasterNodeIterator(typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListItr,
                     typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListBegin,
                     typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListEnd,
                     std::vector<int>::const_iterator IDItr);
#endif
  MasterNodeIterator(const MasterNodeIterator& rhs);

  // Destructor.
  virtual ~MasterNodeIterator();

  // Increment the iterator.
  MasterNodeIterator& operator++();

  // Valid test.
  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Internal method to initialize the state.
  void initialize(typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListItr,
                  typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListBegin,
                  typename std::vector<NodeSpace::NodeList<Dimension>*>::const_iterator nodeListEnd,
                  std::vector<int>::const_iterator IDItr);

#ifndef __GCCXML__
  // The iterator to the current master node ID.
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
#include "MasterNodeIteratorInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class MasterNodeIterator;
}

#endif
