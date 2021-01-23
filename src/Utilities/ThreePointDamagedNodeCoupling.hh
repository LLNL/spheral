//---------------------------------Spheral++----------------------------------//
// ThreePointDamagedNodeCoupling
//
// A functor class encapsulating how we couple solid nodes in the presence of
// multiple materials and damage.
// 
// This for uses the "three point" formalism, which allows damaged points to
// cut communication between pairs that talk across them.
//
// Unlike other NodeCouping objects, ThreePointDamagedNodeCoupling directly
// sets the pairwise f_couple in the NodePairList, so all work is done in the
// constructor.
//
// Created by JMO, Fri Jan  1 15:18:38 PST 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_ThreePointDamagedNodeCoupling__
#define __Spheral_ThreePointDamagedNodeCoupling__

#include "Utilities/NodeCoupling.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

template<typename Dimension>
class ThreePointDamagedNodeCoupling: public NodeCoupling {
public:
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor ;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor.
  ThreePointDamagedNodeCoupling(const FieldList<Dimension, Vector>& position,
                                const FieldList<Dimension, SymTensor>& H,
                                const FieldList<Dimension, SymTensor>& damage,
                                const TableKernel<Dimension>& W,
                                const ConnectivityMap<Dimension>& connectivity,
                                const bool useIntersectConnectivity,
                                NodePairList& pairs);

  // The coupling operator.
  virtual double operator()(const NodePairIdxType& pair) const override {
    return pair.f_couple;
  }

};

}

#endif
