//---------------------------------Spheral++----------------------------------//
// SVPHMassDensityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the mass density in the state according to the SVPH 
// formalism.
//
// Created by JMO, Sat Aug 10 18:52:20 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SVPHMassDensityPolicy_hh__
#define __Spheral_SVPHMassDensityPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class SVPHMassDensityPolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  SVPHMassDensityPolicy(const Scalar& rhoMin,
                        const Scalar& rhoMax);
  virtual ~SVPHMassDensityPolicy();
  
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
  Scalar mRhoMin, mRhoMax;

  SVPHMassDensityPolicy();
  SVPHMassDensityPolicy(const SVPHMassDensityPolicy& rhs);
  SVPHMassDensityPolicy& operator=(const SVPHMassDensityPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SVPHMassDensityPolicy;
}

#endif
