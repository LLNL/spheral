#ifndef MHD_CurrentDensityUPDATEPOLICY_HH
#define MHD_CurrentDensityUPDATEPOLICY_HH

#include <string>
#include "DataBase/UpdatePolicyBase.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension> class DataBase;

//! \class CurrentDensityUpdatePolicy
//! This class alters the sound speed to account for MHD waves.  Currently,
//! it alters the sound speed to reflect the fast magnetosonic wave, which
//! is the faster hydromagnetic propagation.  The slow magnetosonic wave is,
//! hence, not represented.
class CurrentDensityUpdatePolicy:
   public Spheral::UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Vector> >
{
   public:

   //! Construct the update policy.
   //! \param kernel The kernel used for computing SPH interpolants.
   //! \param dataBase The database containing the connectivity data.
   //! \param mu0 The permeability of free space in vacuum.
   CurrentDensityUpdatePolicy(const TableKernel<Dim<3> >& kernel,
                              const DataBase<Dim<3> >& dataBase,
                              double mu0);

   //! Destructor.
   virtual ~CurrentDensityUpdatePolicy();

   //! Overload the methods describing how to update Fields.
   virtual void update(const KeyType& key, 
                       State<Dim<3> >& state,
                       StateDerivatives<Dim<3> >& derivs, 
                       const double multiplier,
                       const double t,
                       const double dt); 

   //! Equivalence.
   virtual bool operator==(const Spheral::UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Vector> >& rhs) const;

   private:

   // Disabled methods.
   CurrentDensityUpdatePolicy(const CurrentDensityUpdatePolicy& rhs);
   CurrentDensityUpdatePolicy& operator=(const CurrentDensityUpdatePolicy& rhs);

   // Kernel used for computing SPH interpolants.
   const TableKernel<Dim<3> >& mKernel;

   // Database used for connectivity data.
   const DataBase<Dim<3> >& mDataBase;

   // Permeability of free space in vacuum.
   double mMu0;
};

}

#else

// Forward declaration.
namespace Spheral {
  class CurrentDensityUpdatePolicy;
}

#endif
