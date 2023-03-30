#ifndef MHD_MAGNETOSONICSPEEDUPDATEPOLICY_HH
#define MHD_MAGNETOSONICSPEEDUPDATEPOLICY_HH

#include <string>
#include "DataBase/UpdatePolicyBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension, typename DataType> class Field;

//! \class MagnetosonicSpeedUpdatePolicy
//! This class alters the sound speed to account for MHD waves.  Currently,
//! it alters the sound speed to reflect the fast magnetosonic wave, which
//! is the faster hydromagnetic propagation.  The slow magnetosonic wave is,
//! hence, not represented.
class MagnetosonicSpeedUpdatePolicy:
   public Spheral::UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >
{
   public:

   //! Construct the update policy.
   //! \param mu0 The permeability of free space in vacuum.
   explicit MagnetosonicSpeedUpdatePolicy(double mu0);

   //! Destructor.
   virtual ~MagnetosonicSpeedUpdatePolicy();

   //! Overload the methods describing how to update Fields.
   virtual void update(const KeyType& key, 
                       State<Dim<3> >& state,
                       StateDerivatives<Dim<3> >& derivs, 
                       const double multiplier,
                       const double t,
                       const double dt); 

   //! Equivalence.
   virtual bool operator==(const Spheral::UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >& rhs) const;

   private:

   // Disabled methods.
   MagnetosonicSpeedUpdatePolicy(const MagnetosonicSpeedUpdatePolicy& rhs);
   MagnetosonicSpeedUpdatePolicy& operator=(const MagnetosonicSpeedUpdatePolicy& rhs);

   // Permeability of free space in vacuum.
   double mMu0;
};

}
}

#else

// Forward declaration.
namespace Spheral {
  class MagnetosonicSpeedUpdatePolicy;
}

#endif
