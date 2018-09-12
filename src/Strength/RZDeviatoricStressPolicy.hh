//---------------------------------Spheral++----------------------------------//
// RZDeviatoricStressPolicy.
//
// Specialized version of the IncrementFieldList policy, with some criteria for
// zeroing out the deviatoric stress in special cases.
// This one is also specialized for the RZ strength case.
//
// Created by JMO, Mon May  9 14:09:12 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_RZDeviatoricStress_hh__
#define __Spheral_RZDeviatoricStress_hh__

#include "DataBase/IncrementFieldList.hh"
#include "Geometry/Dimension.hh"

#include <string>

namespace Spheral {

class RZDeviatoricStressPolicy: public IncrementFieldList<Dim<2>, Dim<2>::SymTensor> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  RZDeviatoricStressPolicy();
  virtual ~RZDeviatoricStressPolicy();
  
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
  RZDeviatoricStressPolicy(const RZDeviatoricStressPolicy& rhs);
  RZDeviatoricStressPolicy& operator=(const RZDeviatoricStressPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  class RZDeviatoricStressPolicy;
}

#endif
