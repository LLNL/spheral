//---------------------------------Spheral++----------------------------------//
// JohnsonCookDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the tensor damage variable with the Johnson-Cook model.
//
// Created by JMO, Thu Jul 12 15:33:26 PDT 2018
//----------------------------------------------------------------------------//
#ifndef __Spheral_JohnsonCookDamagePolicy_hh__
#define __Spheral_JohnsonCookDamagePolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension> class JohnsonCookDamageModel;

template<typename Dimension>
class JohnsonCookDamagePolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  explicit JohnsonCookDamagePolicy();
  virtual ~JohnsonCookDamagePolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  JohnsonCookDamagePolicy(const JohnsonCookDamagePolicy& rhs);
  JohnsonCookDamagePolicy& operator=(const JohnsonCookDamagePolicy& rhs);
};

}

#endif
