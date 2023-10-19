//---------------------------------Spheral++----------------------------------//
// PorousBulkModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent bulk modulus state in the presence of porosity.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorousBulkModulusPolicy_hh__
#define __Spheral_PorousBulkModulusPolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class PorousBulkModulusPolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;

  // Constructors, destructor.
  PorousBulkModulusPolicy();
  virtual ~PorousBulkModulusPolicy();
  
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
  PorousBulkModulusPolicy(const PorousBulkModulusPolicy& rhs);
  PorousBulkModulusPolicy& operator=(const PorousBulkModulusPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PorousBulkModulusPolicy;
}

#endif
