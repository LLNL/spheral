//---------------------------------Spheral++----------------------------------//
// MeshIdealHPolicy
//
// Created by JMO, Fri Aug 16 14:47:57 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_MeshIdealHPolicy_hh__
#define __Spheral_MeshIdealHPolicy_hh__

#include <string>

#include "DataBase/ReplaceBoundedState.hh"
#include "NodeList/SmoothingScaleBase.hh"

namespace Spheral {

template<typename Dimension>
class MeshIdealHPolicy: 
    public ReplaceBoundedState<Dimension, typename Dimension::SymTensor, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename FieldUpdatePolicyBase<Dimension, SymTensor>::KeyType KeyType;

  // Constructors, destructor.
  MeshIdealHPolicy(const SmoothingScaleBase<Dimension>& smoothingScaleBase,
                   const Scalar hmin,
                   const Scalar hmax,
                   const Scalar hminratio,
                   const Scalar nPerh);
  virtual ~MeshIdealHPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // If the derivative stored values for the pair-accelerations has not been updated,
  // we need to just time advance normally.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) {
    ReplaceBoundedState<Dimension, SymTensor, Scalar>::update(key,
                                                              state,
                                                              derivs,
                                                              multiplier,
                                                              t,
                                                              dt);
  }

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  static bool mFired;
  const SmoothingScaleBase<Dimension>& mSmoothingScaleBase;
  Scalar mhmin, mhmax, mhminratio, mnPerh;

  MeshIdealHPolicy(const MeshIdealHPolicy& rhs);
  MeshIdealHPolicy& operator=(const MeshIdealHPolicy& rhs);
};

}

#endif
