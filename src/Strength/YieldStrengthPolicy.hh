//---------------------------------Spheral++----------------------------------//
// YieldStrengthPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent yield strength state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_YieldStrengthPolicy_hh__
#define __Spheral_YieldStrengthPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class YieldStrengthPolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  YieldStrengthPolicy(const bool scaleWithPorosity = false);
  virtual ~YieldStrengthPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  // Should the yield strength be scaled with porosity?
  bool scaleWithPorosity() const       { return mScaleWithPorosity; }
  void scaleWithPorosity(const bool x) { mScaleWithPorosity = x; }

private:
  //--------------------------- Private Interface ---------------------------//
  bool mScaleWithPorosity;

  YieldStrengthPolicy(const YieldStrengthPolicy& rhs);
  YieldStrengthPolicy& operator=(const YieldStrengthPolicy& rhs);
};

}

#endif
