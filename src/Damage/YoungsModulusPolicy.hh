//---------------------------------Spheral++----------------------------------//
// YoungsModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent Youngs modulus.
//
// Created by JMO, Sun Jan 23 14:14:24 PST 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_YoungsModulusPolicy_hh__
#define __Spheral_YoungsModulusPolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class YoungsModulusPolicy: 
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
  YoungsModulusPolicy();
  virtual ~YoungsModulusPolicy();
  
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
  YoungsModulusPolicy(const YoungsModulusPolicy& rhs);
  YoungsModulusPolicy& operator=(const YoungsModulusPolicy& rhs);
};

}

#endif
