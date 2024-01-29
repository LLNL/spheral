//---------------------------------Spheral++----------------------------------//
// HullVolumePolicy -- An implementation of ReplaceState specialized
// for the updating the volume based on the local hull constructions.
//
// Created by JMO, Wed Jul  2 22:35:26 PDT 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_HullVolumePolicy_hh__
#define __Spheral_HullVolumePolicy_hh__

#include <string>

#include "DataBase/ReplaceFieldList.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

template<typename Dimension>
class HullVolumePolicy: public ReplaceFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  HullVolumePolicy(const ConnectivityMap<Dimension>& connectivityMap);
  virtual ~HullVolumePolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // We'll make the updateAsIncrement a no-op.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) {}

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  const ConnectivityMap<Dimension>& mConnectivityMap;

  HullVolumePolicy(const HullVolumePolicy& rhs);
  HullVolumePolicy& operator=(const HullVolumePolicy& rhs);
};

}

#endif
