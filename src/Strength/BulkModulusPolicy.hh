//---------------------------------Spheral++----------------------------------//
// BulkModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent bulk modulus state.
//
// Created by JMO, Tue Oct 5 16:40:54 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_BulkModulusPolicy_hh__
#define __Spheral_BulkModulusPolicy_hh__

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

template<typename Dimension>
class BulkModulusPolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  BulkModulusPolicy();
  virtual ~BulkModulusPolicy();
  
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
  BulkModulusPolicy(const BulkModulusPolicy& rhs);
  BulkModulusPolicy& operator=(const BulkModulusPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class BulkModulusPolicy;
}

#endif
