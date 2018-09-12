//---------------------------------Spheral++----------------------------------//
// RZPlasticStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent plastic strain state.
// This one is also specialized for the RZ strength case.
//
// Created by JMO, Mon May  9 14:09:12 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_RZPlasticStrainPolicy_hh__
#define __Spheral_RZPlasticStrainPolicy_hh__

#include "DataBase/FieldListUpdatePolicyBase.hh"
#include "Geometry/Dimension.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

class RZPlasticStrainPolicy: 
    public FieldListUpdatePolicyBase<Dim<2>, Dim<2>::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  RZPlasticStrainPolicy();
  virtual ~RZPlasticStrainPolicy();
  
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
  RZPlasticStrainPolicy(const RZPlasticStrainPolicy& rhs);
  RZPlasticStrainPolicy& operator=(const RZPlasticStrainPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  class RZPlasticStrainPolicy;
}

#endif
