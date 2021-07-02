//---------------------------------Spheral++----------------------------------//
// PairMaxDamageNodeCoupling
//
// A functor class encapsulating how we couple solid nodes in the presence of
// multiple materials and damage.
//
// This form simply directly damages each pair based on their mutual damage.
//
// Created by JMO: Fri Jul 31 14:46:25 PDT 2015
//        updated: Fri Feb  5 12:48:20 PST 2021 JMO
//----------------------------------------------------------------------------//
#ifndef __Spheral_PairMaxDamageNodeCoupling__
#define __Spheral_PairMaxDamageNodeCoupling__

#include "Utilities/NodeCoupling.hh"
#include "DataBase/State.hh"

namespace Spheral {

template<typename Dimension>
class PairMaxDamageNodeCoupling: public NodeCoupling {
public:
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructor.
  PairMaxDamageNodeCoupling(const State<Dimension>& state,
                      NodePairList& pairs);

private:
  // Forbidden methods.
  PairMaxDamageNodeCoupling();
};

}

#endif
