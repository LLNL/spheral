//---------------------------------Spheral++----------------------------------//
// DamagedNodeCoupling
//
// A functor class encapsulating how we couple solid nodes in the presence of
// multiple materials and damage.
//
// Created by JMO, Fri Jul 31 14:46:25 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_DamagedNodeCoupling__
#define __Spheral_DamagedNodeCoupling__

#include "Utilities/NodeCoupling.hh"
#include "Utilities/DBC.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {

template<typename Dimension>
class DamagedNodeCoupling: public NodeCoupling {
public:
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructor.
  DamagedNodeCoupling(const FieldList<Dimension, SymTensor>& damage,
                      NodePairList& pairs):
    NodeCoupling(),
    mDamage(damage) {
    const auto n = pairs.size();
#pragma omp for
    for (auto k = 0u; k < n; ++k) {
      pairs[k].f_couple = (*this)(pairs[k]);
    }
  }

  // The coupling operator.
  virtual double operator()(const NodePairIdxType& pair) const override {
    const auto sDi = scalarDamage(pair.i_list, pair.i_node);
    const auto sDj = scalarDamage(pair.j_list, pair.j_node);
    return 1.0 - std::max(sDi, sDj);
  }

private:
  const FieldList<Dimension, SymTensor>& mDamage;

  // Extract the effective scalar damage on a node.
  double scalarDamage(const unsigned nodeListi, const unsigned i) const {
    Scalar sDi = std::max(0.0, std::min(1.0, mDamage(nodeListi, i).eigenValues().maxElement()));
    if (sDi > 1.0 - 1.0e-3) sDi = 1.0;
    return sDi;
  }

  // // Construct a unit vector of the argument, going to zero as the magnitude falls below a given "fuzz".
  // Vector unitVectorWithZero(const Vector& x,
  //                           const double fuzz = 0.01) const {
  //   if (x.magnitude2() < fuzz) {
  //     return Vector::zero;
  //   } else {
  //     return x.unitVector();
  //   }
  // }

  // // Damage direction based on the gradient.
  // Vector damageDirection(const unsigned nodeListi, const unsigned i) const {
  //   return unitVectorWithZero(mDamageGradient(nodeListi, i)*Dimension::nDim/mH(nodeListi, i).Trace());
  // }

  // Forbidden methods.
  DamagedNodeCoupling();
};

}

#endif
