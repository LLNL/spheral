//---------------------------------Spheral++----------------------------------//
// ALEPositionPolicy -- This is basically a direct copy of the standard 
//                      position policy but instead we're substituting in 
//                      the nodal velocity as the derivative.
//
// J. M. Pearl 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_ALEPositionPolicy_hh__
#define __Spheral_ALEPositionPolicy_hh__

#include "DataBase/IncrementFieldList.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class ALEPositionPolicy: 
    public IncrementFieldList<Dimension, typename Dimension::Vector> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename FieldListUpdatePolicyBase<Dimension, Vector>::KeyType KeyType;

  // Constructors, destructor.
  ALEPositionPolicy();
  virtual ~ALEPositionPolicy();
  
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
  ALEPositionPolicy(const ALEPositionPolicy& rhs);
  ALEPositionPolicy& operator=(const ALEPositionPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ALEPositionPolicy;
}

#endif
