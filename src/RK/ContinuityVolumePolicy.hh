//---------------------------------Spheral++----------------------------------//
// ContinuityVolumePolicy -- An implementation of IncrementFieldList
// specialized for time evolving the volume per point using the continuity
// equation.
//
// Created by JMO, Tue Sep 20 14:53:32 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_ContinuityVolumePolicy_hh__
#define __Spheral_ContinuityVolumePolicy_hh__

#include <string>

#include "DataBase/IncrementFieldList.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension>
class ContinuityVolumePolicy: public IncrementFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  ContinuityVolumePolicy();
  virtual ~ContinuityVolumePolicy();
  
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
  ContinuityVolumePolicy(const ContinuityVolumePolicy& rhs);
  ContinuityVolumePolicy& operator=(const ContinuityVolumePolicy& rhs);
};

}

#endif
