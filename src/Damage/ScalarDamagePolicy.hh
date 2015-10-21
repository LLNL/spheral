//---------------------------------Spheral++----------------------------------//
// ScalarDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent scalar damage state.
//
// Created by JMO, Sun Sep 26 16:15:17 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ScalarDamagePolicy_hh__
#define __Spheral_ScalarDamagePolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
namespace NodeSpace {
  template<typename Dimension> class FluidNodeList;
}
namespace FieldSpace {
  template<typename Dimension, typename DataType> class Field;
}
namespace PhysicsSpace {
  template<typename Dimension> class DamageModel;
}

template<typename Dimension>
class ScalarDamagePolicy: 
    public UpdatePolicyBase<Dimension, FieldSpace::Field<Dimension, typename Dimension::Scalar> > {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename FieldSpace::Field<Dimension, Scalar> FieldType;
  typedef typename UpdatePolicyBase<Dimension, FieldType>::KeyType KeyType;

  // Constructors, destructor.
  explicit ScalarDamagePolicy(const PhysicsSpace::DamageModel<Dimension>& damageModel);
  virtual ~ScalarDamagePolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension, FieldType>& rhs) const;

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  const PhysicsSpace::DamageModel<Dimension>* mDamageModelPtr;

  ScalarDamagePolicy();
  ScalarDamagePolicy(const ScalarDamagePolicy& rhs);
  ScalarDamagePolicy& operator=(const ScalarDamagePolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ScalarDamagePolicy;
}

#endif
