//---------------------------------Spheral++----------------------------------//
// NodeCoupling
//
// A functor base class encapsulating how we couple pairs of nodes.
//
// Created by JMO, Fri Jul 31 21:16:59 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeCoupling__
#define __Spheral_NodeCoupling__

#include "Neighbor/NodePairList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// The generic node coupling -- all nodes coupled equally.
class NodeCoupling {
public:
  // Constructor, destructor.
  NodeCoupling() {}
  virtual ~NodeCoupling() {}
  

  // The coupling operator.
  virtual double operator()(const NodePairIdxType& /*pair*/) const {
    return 1.0;
  }
};

// A variant where only nodes within a NodeList are coupled.
class PerNodeListNodeCoupling : public NodeCoupling {
public:
  // Constructor.
  PerNodeListNodeCoupling(): NodeCoupling() {}
  virtual ~PerNodeListNodeCoupling() {}

  // The coupling operator.
  virtual double operator()(const NodePairIdxType& pair) const override {
    return (pair.i_list == pair.j_list ? 1.0 : 0.0);
  }
};

}

#endif
