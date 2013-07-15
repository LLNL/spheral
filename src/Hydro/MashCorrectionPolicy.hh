//---------------------------------Spheral++----------------------------------//
// MashCorrectionPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the MASH node based correction.
//
// Created by JMO, Mon Sep 26 21:23:07 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_MashCorrectionPolicy_hh__
#define __Spheral_MashCorrectionPolicy_hh__

#include <string>

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
class MashCorrectionPolicy: 
    public UpdatePolicyBase<Dimension, FieldSpace::Field<Dimension, typename Dimension::Tensor> > {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Tensor Tensor;
  typedef typename FieldSpace::Field<Dimension, Tensor> FieldType;
  typedef typename UpdatePolicyBase<Dimension, FieldType>::KeyType KeyType;

  // Constructors, destructor.
  MashCorrectionPolicy();
  virtual ~MashCorrectionPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension, FieldType>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  MashCorrectionPolicy(const MashCorrectionPolicy& rhs);
  MashCorrectionPolicy& operator=(const MashCorrectionPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class MashCorrectionPolicy;
}

#endif
