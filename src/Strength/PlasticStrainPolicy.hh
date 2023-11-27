//---------------------------------Spheral++----------------------------------//
// PlasticStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent plastic strain state.
//
// Created by JMO, Wed Sep 15 23:07:21 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PlasticStrainPolicy_hh__
#define __Spheral_PlasticStrainPolicy_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class PlasticStrainPolicy: 
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
  PlasticStrainPolicy();
  virtual ~PlasticStrainPolicy();
  
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
  PlasticStrainPolicy(const PlasticStrainPolicy& rhs);
  PlasticStrainPolicy& operator=(const PlasticStrainPolicy& rhs);
};

}

#endif
