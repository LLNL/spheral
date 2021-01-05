//---------------------------------Spheral++----------------------------------//
// DamagedNodeCouplingWithFrags
//
// A functor class encapsulating how we couple solid nodes in the presence of
// multiple materials and damage.  This version adds logic to decouple based
// on fragment ID as well.
//
// Created by JMO, Fri Aug  7 09:04:16 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_DamagedNodeCouplingWithFrags__
#define __Spheral_DamagedNodeCouplingWithFrags__

#include "Utilities/DamagedNodeCoupling.hh"
#include "Utilities/DBC.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {

template<typename Dimension>
class DamagedNodeCouplingWithFrags: public DamagedNodeCoupling<Dimension> {
public:
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructor.
  DamagedNodeCouplingWithFrags(const FieldList<Dimension, SymTensor>& damage,
                               const FieldList<Dimension, Vector>& damageGradient,
                               const FieldList<Dimension, SymTensor>& H,
                               const FieldList<Dimension, int>& fragIDs):
    DamagedNodeCoupling<Dimension>(damage, damageGradient, H),
    mFragIDs(fragIDs) {}

  // The coupling operator.
  virtual double operator()(const NodePairIdxType& pair) const override {
    if (mFragIDs(pair.i_list, pair.i_node) != mFragIDs(pair.j_list, pair.j_node)) {
      return 0.0;
    } else {
      return DamagedNodeCoupling<Dimension>::operator()(pair);
    }
  }

private:
  const FieldList<Dimension, int>& mFragIDs;

  // Forbidden methods.
  DamagedNodeCouplingWithFrags();
};

}

#endif
