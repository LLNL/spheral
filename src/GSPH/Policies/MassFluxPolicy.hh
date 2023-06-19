//---------------------------------Spheral++----------------------------------//
// MassFluxPolicy -- This is basically a direct copy of the standard 
//                      position policy but instead we're substituting in 
//                      the nodal velocity as the derivative.
//
// J. M. Pearl 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_MassFluxPolicy_hh__
#define __Spheral_MassFluxPolicy_hh__

#include "DataBase/IncrementFieldList.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class MassFluxPolicy: 
    public IncrementFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  MassFluxPolicy();
  MassFluxPolicy(const std::string& depend0);
  MassFluxPolicy(const std::string& depend0,const std::string& depend1);
  MassFluxPolicy(const std::string& depend0,const std::string& depend1,const std::string& depend2);
  MassFluxPolicy(const std::string& depend0,const std::string& depend1,const std::string& depend2,const std::string& depend3);
  MassFluxPolicy(const std::string& depend0,const std::string& depend1,const std::string& depend2,const std::string& depend3,const std::string& depend4);
  MassFluxPolicy(const std::string& depend0,const std::string& depend1,const std::string& depend2,const std::string& depend3,const std::string& depend4,const std::string& depend5);
  virtual ~MassFluxPolicy();
  
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
  MassFluxPolicy(const MassFluxPolicy& rhs);
  MassFluxPolicy& operator=(const MassFluxPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class MassFluxPolicy;
}

#endif
