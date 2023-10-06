//---------------------------------Spheral++----------------------------------//
// ContinuityVolumePolicy -- An implementation of IncrementFieldList
// specialized for time evolving the volume per point using the continuity
// equation.
//
// Specialized for RZ geometry.
//
// Created by JMO, Tue Oct 29 16:14:03 PDT 2019
//----------------------------------------------------------------------------//
#ifndef __Spheral_ContinuityVolumePolicyRZ_hh__
#define __Spheral_ContinuityVolumePolicyRZ_hh__

#include <string>

#include "DataBase/IncrementState.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class ContinuityVolumePolicyRZ: public IncrementState<Dim<2>, Dim<2>::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using SymTensor = Dimension::SymTensor;
  using KeyType = UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  ContinuityVolumePolicyRZ();
  virtual ~ContinuityVolumePolicyRZ();
  
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
  ContinuityVolumePolicyRZ(const ContinuityVolumePolicyRZ& rhs);
  ContinuityVolumePolicyRZ& operator=(const ContinuityVolumePolicyRZ& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  class ContinuityVolumePolicyRZ;
}

#endif
