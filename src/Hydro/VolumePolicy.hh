//---------------------------------Spheral++----------------------------------//
// VolumePolicy -- An implementation of ReplaceState specialized
// for the updating the volume based on the current Voronoi tesselation.
//
// Created by JMO, Mon Aug  1 15:17:36 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_VolumePolicy_hh__
#define __Spheral_VolumePolicy_hh__

#include <string>

#include "DataBase/ReplaceFieldList.hh"

namespace Spheral {

template<typename Dimension>
class VolumePolicy: public ReplaceFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  VolumePolicy();
  virtual ~VolumePolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // We'll make the updateAsIncrement a no-op.
  virtual void updateAsIncrement(const KeyType& /*key*/,
                                 State<Dimension>& /*state*/,
                                 StateDerivatives<Dimension>& /*derivs*/,
                                 const double /*multiplier*/,
                                 const double /*t*/,
                                 const double /*dt*/) {}

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  VolumePolicy(const VolumePolicy& rhs);
  VolumePolicy& operator=(const VolumePolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class VolumePolicy;
}

#endif
