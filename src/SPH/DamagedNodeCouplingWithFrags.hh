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

#include "DamagedNodeCoupling.hh"
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
  DamagedNodeCouplingWithFrags(const FieldSpace::FieldList<Dimension, SymTensor>& damage,
                               const FieldSpace::FieldList<Dimension, Vector>& damageGradient,
                               const FieldSpace::FieldList<Dimension, SymTensor>& H,
                               const FieldSpace::FieldList<Dimension, int>& fragIDs):
    DamagedNodeCoupling<Dimension>(damage, damageGradient, H),
    mFragIDs(fragIDs) {}

  // The coupling operator.
  virtual double operator()(const unsigned nodeListi, const unsigned i,
                            const unsigned nodeListj, const unsigned j) const {
    if (mFragIDs(nodeListi, i) != mFragIDs(nodeListj, j)) {
      return 0.0;
    } else {
      return DamagedNodeCoupling<Dimension>::operator()(nodeListi, i, nodeListj, j);
    }
  }

private:
  const FieldSpace::FieldList<Dimension, int>& mFragIDs;

  // Forbidden methods.
  DamagedNodeCouplingWithFrags();
};

}

#endif
