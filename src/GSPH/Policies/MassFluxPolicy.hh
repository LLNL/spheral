//---------------------------------Spheral++----------------------------------//
// MassFluxPolicy -- update method for ALE - based hydro schemes that allow
//                   for mass flux between nodes.
//
// J. M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_MassFluxPolicy_hh__
#define __Spheral_MassFluxPolicy_hh__

#include "DataBase/IncrementState.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class MassFluxPolicy: 
    public IncrementState<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename IncrementState<Dimension,Scalar>::KeyType;

  // Constructors, destructor.
  MassFluxPolicy(std::initializer_list<std::string> depends = {});
  virtual ~MassFluxPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const  override;

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
