//---------------------------------Spheral++----------------------------------//
// JohnsonCookFailureStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the strain.
//
// Created by JMO, Thu Jul 12 13:40:45 PDT 2018
//----------------------------------------------------------------------------//
#ifndef __Spheral_JohnsonCookFailureStrainPolicy_hh__
#define __Spheral_JohnsonCookFailureStrainPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension> class StrainModel;

template<typename Dimension>
class JohnsonCookFailureStrainPolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef Field<Dimension, Scalar> FieldType;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  JohnsonCookFailureStrainPolicy(const Field<Dimension, Scalar>& D1,
                                 const Field<Dimension, Scalar>& D2,
                                 const double D3,
                                 const double D4,
                                 const double D5,
                                 const double epsilondot0,
                                 const double sigmamax,
                                 const double efailmin,
                                 const double Tcrit);
  virtual ~JohnsonCookFailureStrainPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  const Field<Dimension, Scalar>& mD1;
  const Field<Dimension, Scalar>& mD2;
  double mD3, mD4, mD5, mepsilondot0, msigmamax, mefailmin, mTcrit;

  // No copy or default construction.
  JohnsonCookFailureStrainPolicy();
  JohnsonCookFailureStrainPolicy(const JohnsonCookFailureStrainPolicy& rhs);
  JohnsonCookFailureStrainPolicy& operator=(const JohnsonCookFailureStrainPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class JohnsonCookFailureStrainPolicy;
}

#endif
