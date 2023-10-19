//---------------------------------Spheral++----------------------------------//
// PorositySolidMassDensityPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the dependent solid mass density in the presence
// of porosity.
//
// Created by JMO, Wed Oct  4 13:51:43 PDT 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorositySolidMassDensityPolicy_hh__
#define __Spheral_PorositySolidMassDensityPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class PorositySolidMassDensityPolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  PorositySolidMassDensityPolicy();
  virtual ~PorositySolidMassDensityPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  PorositySolidMassDensityPolicy(const PorositySolidMassDensityPolicy& rhs);
  PorositySolidMassDensityPolicy& operator=(const PorositySolidMassDensityPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PorositySolidMassDensityPolicy;
}

#endif
