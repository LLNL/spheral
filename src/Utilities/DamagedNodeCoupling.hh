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
                      const FieldList<Dimension, Vector>& damageGradient,
                      const FieldList<Dimension, SymTensor>& H):
    NodeCoupling(),
    mDamage(damage),
    mDamageGradient(damageGradient),
    mH(H) {}

  // The coupling operator.
  virtual double operator()(const unsigned nodeListi, const unsigned i,
                            const unsigned nodeListj, const unsigned j) const {
    const Scalar sDi = scalarDamage(nodeListi, i);
    const Scalar sDj = scalarDamage(nodeListj, j);
    if (sDi + sDj > 1.95) {
      return 1.0 - std::max(sDi, sDj);
    } else {
      const Vector gradDi = damageDirection(nodeListi, i);
      const Vector gradDj = damageDirection(nodeListj, j);
      const Scalar gradDdot = gradDi.dot(gradDj);
      const Scalar phi = ((std::abs(gradDdot) < 0.1 or std::min(sDi, sDj) < 0.05) ? 
                          1.0 :
                          std::max(0.0, std::min(1.0, -gradDdot)));
      CHECK(phi >= 0.0 and phi <= 1.0);
      const Scalar fDeffij = FastMath::pow4(std::max(0.0, std::min(1.0, 1.0 - phi*std::max(sDi, sDj))));
      ENSURE(fDeffij >= 0.0 and fDeffij <= 1.0);
      return fDeffij;
    }
  }

private:
  const FieldList<Dimension, SymTensor>& mDamage;
  const FieldList<Dimension, Vector>& mDamageGradient;
  const FieldList<Dimension, SymTensor>& mH;

  // Extract the effective scalar damage on a node.
  double scalarDamage(const unsigned nodeListi, const unsigned i) const {
    Scalar sDi = std::max(0.0, std::min(1.0, mDamage(nodeListi, i).eigenValues().maxElement()));
    if (sDi > 1.0 - 1.0e-3) sDi = 1.0;
    return sDi;
  }

  // Construct a unit vector of the argument, going to zero as the magnitude falls below a given "fuzz".
  Vector unitVectorWithZero(const Vector& x,
                            const double fuzz = 0.01) const {
    if (x.magnitude2() < fuzz) {
      return Vector::zero;
    } else {
      return x.unitVector();
    }
  }

  // Damage direction based on the gradient.
  Vector damageDirection(const unsigned nodeListi, const unsigned i) const {
    return unitVectorWithZero(mDamageGradient(nodeListi, i)*Dimension::nDim/mH(nodeListi, i).Trace());
  }

  // Forbidden methods.
  DamagedNodeCoupling();
};

}

#endif
