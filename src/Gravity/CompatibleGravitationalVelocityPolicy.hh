//---------------------------------Spheral++----------------------------------//
// CompatibleGravitationalVelocityPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the gravitationally driven velocity in an exactly energy conservative manner.
// 
// Created by JMO, Mon Oct  2 13:44:59 PDT 2017
//----------------------------------------------------------------------------//

#ifndef __Spheral_CompatibleGravitationalVelocityPolicy_hh__
#define __Spheral_CompatibleGravitationalVelocityPolicy_hh__

#include <string>

#include "DataBase/IncrementFieldList.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
namespace NodeSpace {
  template<typename Dimension> class FluidNodeList;
}
namespace FieldSpace {
  template<typename Dimension, typename DataType> class FieldList;
}
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}

template<typename Dimension>
class CompatibleGravitationalVelocityPolicy: 
    public IncrementFieldList<Dimension, typename Dimension::Vector> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename FieldListUpdatePolicyBase<Dimension, Vector>::KeyType KeyType;

  // Constructors, destructor.
  CompatibleGravitationalVelocityPolicy(const DataBaseSpace::DataBase<Dimension>& db,
                                        const Scalar G,
                                        const Scalar softeningLength);
  virtual ~CompatibleGravitationalVelocityPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // For intra-timestep velocity estimates, just use the regular time advancement.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) {
    IncrementFieldList<Dimension, Vector>::update(key,
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
  const DataBaseSpace::DataBase<Dimension>* mDataBasePtr;
  const Scalar mG, mSofteningLength;
  FieldSpace::FieldList<Dimension, Vector> mPositions0, mVelocity0;

  CompatibleGravitationalVelocityPolicy(const CompatibleGravitationalVelocityPolicy& rhs);
  CompatibleGravitationalVelocityPolicy& operator=(const CompatibleGravitationalVelocityPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class CompatibleGravitationalVelocityPolicy;
}

#endif
