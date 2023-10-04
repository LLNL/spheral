//---------------------------------Spheral++----------------------------------//
// PalphaPressurePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent pressure state for use with the P-alpha
// porosity model.
//
// Created by JMO, Wed Aug 30 13:36:37 PDT 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_PalphaPressurePolicy_hh__
#define __Spheral_PalphaPressurePolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;

template<typename Dimension>
class PalphaPressurePolicy: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename FieldUpdatePolicy<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  PalphaPressurePolicy();
  virtual ~PalphaPressurePolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  PalphaPressurePolicy(const PalphaPressurePolicy& rhs);
  PalphaPressurePolicy& operator=(const PalphaPressurePolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PalphaPressurePolicy;
}

#endif
