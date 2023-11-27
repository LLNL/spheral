//---------------------------------Spheral++----------------------------------//
// MashCorrectionPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the MASH node based correction.
//
// Created by JMO, Mon Sep 26 21:23:07 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_MashCorrectionPolicy_hh__
#define __Spheral_MashCorrectionPolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class MashCorrectionPolicy: 
    public UpdatePolicyBase<Dimension, Field<Dimension, typename Dimension::Tensor> > {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Tensor Tensor;
  typedef typename Field<Dimension, Tensor> FieldType;
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

#endif
