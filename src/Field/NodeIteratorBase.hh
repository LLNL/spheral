//---------------------------------Spheral++----------------------------------//
// NodeIteratorBase -- This is a helper class for the DataBase.  This is the
// base class for the node iterators.  This class implements the basic methods,
// most importantly how to index by such iterators into Fields/FieldLists.
// Based on the old NodeIDIterator class.
//
// Created by J. Michael Owen, Mon Mar 17 13:11:19 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeIteratorBase_hh__
#define __Spheral_NodeIteratorBase_hh__

#include <vector>

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension> class FluidNodeList;

template<typename Dimension>
class NodeIteratorBase {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors and destructors.
  NodeIteratorBase();
  NodeIteratorBase(const NodeIteratorBase& rhs);

  // Destructor.
  virtual ~NodeIteratorBase();

  // Assignment.
  NodeIteratorBase& operator=(const NodeIteratorBase& rhs);

  // Access the NodeList iterator/pointer.
  typename std::vector<NodeList<Dimension>*>::const_iterator nodeListIterator() const;
  const NodeList<Dimension>* nodeListPtr() const;
  const FluidNodeList<Dimension>* fluidNodeListPtr() const;

  // Access the node ID.
  int nodeID() const;

  // Calculate the index of the Field currently pointed at by the iterator.
  int fieldID() const;

  // Stuff we can extract about the node.
  Scalar mass() const;
  const Vector& position() const;
  const Vector& velocity() const;
  const SymTensor& H() const;

  // Test if the node is internal or a ghost.
  // Note that if this isn't a valid node (such as an end() iterator) it is
  // neither internal nor a ghost.
  bool internalNode() const;
  bool ghostNode() const;

  // Comparison operators.
  bool operator==(const NodeIteratorBase& rhs) const;
  bool operator!=(const NodeIteratorBase& rhs) const;
  bool operator<(const NodeIteratorBase& rhs) const;
  bool operator>(const NodeIteratorBase& rhs) const;
  bool operator<=(const NodeIteratorBase& rhs) const;
  bool operator>=(const NodeIteratorBase& rhs) const;

  // Minimal method to test that the NodeIteratorBase is in a valid state.
  virtual bool valid() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Internal state data.
  int mNodeID;
  int mFieldID;
  typename std::vector<NodeList<Dimension>*>::const_iterator mNodeListBegin;
  typename std::vector<NodeList<Dimension>*>::const_iterator mNodeListEnd;
  typename std::vector<NodeList<Dimension>*>::const_iterator mNodeListItr;
};

}

#include "NodeIteratorBaseInline.hh"

#endif
