//---------------------------------Spheral++----------------------------------//
// PairwisePlasticStrainPolicy -- Modification to the default plastic strain  
//                                policy to allow for some pairwise behavior.
//                                This version stores 1/sqrt(J2) in a state
//                                variable attached to the hydro. The pairwise
//                                stress deviators in the eval derivs loop are 
//                                scaled by the plastic yield factor based on
//                                the pairwise minimum yield strength and their
//                                respective stored J2 values. At damage fronts
//                                the yield is effective the minimum pairwise
//                                yield instead of 1/2 the maximum.
//
//  J M Pearl 08/2023 
//  (essentially copy-paste from JMO's Strength/PlasticStrainPolicy)
//----------------------------------------------------------------------------//
#ifndef __Spheral_PairwisePlasticStrainPolicy_hh__
#define __Spheral_PairwisePlasticStrainPolicy_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class PairwisePlasticStrainPolicy: 
    public FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  PairwisePlasticStrainPolicy();
  virtual ~PairwisePlasticStrainPolicy();
  
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
  PairwisePlasticStrainPolicy(const PairwisePlasticStrainPolicy& rhs);
  PairwisePlasticStrainPolicy& operator=(const PairwisePlasticStrainPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PairwisePlasticStrainPolicy;
}

#endif
