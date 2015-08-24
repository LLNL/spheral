#ifndef MHD_SpitzerResistivityUpdatePolicy_HH
#define MHD_SpitzerResistivityUpdatePolicy_HH

#include <string>
#include "DataBase/UpdatePolicyBase.hh"
#include "Field/Field.hh"
#include "Geometry/Dimension.hh"


namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
namespace FieldSpace {
   template<typename Dimension, typename DataType> class Field;
}

namespace MHDSpace {

//! \class SpitzerResistivityUpdatePolicy
//! The classical Spitzer model for the resistivity of a plasma is 
//!
//!                C
//!       R =    -----,   T  > 0  and  R = Rmax, T  == 0
//!                3/2     e                      e
//!               T
//!                e
//!
//! where C is some given constant, Rmax is the maximum resistivity allowed, 
//! and T  is the electron temperature of the plasma.
//!      e
//!
//! To allow for a simple "anomalous resistivity" model, we can also make the 
//! model return the minimum Resistivity when the mass density falls below 
//! a minimum threshhold, rhoMin.
//! So this model has two parameters: C, Rmax, and rhoMin.
class SpitzerResistivityUpdatePolicy:
   public Spheral::UpdatePolicyBase<Dim<3>, FieldSpace::Field<Dim<3>, Dim<3>::Scalar> >
{
   public:

   //! Construct the Spitzer policy with the given parameters.  This 
   //! constructor sets the minimum density threshhold to 0, effectively 
   //! ignoring "anomalous" Resistivity.
   //! \param C The model constant.
   //! \param Rmax The maximum (parallel) resistivity this model is allowed to return.
   SpitzerResistivityUpdatePolicy(double C, double Rmax = FLT_MAX);

   //! Construct the Spitzer policy with the given parameters AND the minimum density.
   //! \param C The model constant.
   //! \param Rmax The maximum (parallel) resistivity this model is allowed to return.
   //! \param rhoMin The minimum density threshhold (0 by default).
   SpitzerResistivityUpdatePolicy(double C, double Rmax, double rhoMin);

   //! Destructor.
   ~SpitzerResistivityUpdatePolicy();

   //! The update method override.
   void update(const KeyType& key, 
               State<Dim<3> >& state,
               StateDerivatives<Dim<3> >& derivs, 
               const double multiplier,
               const double t,
               const double dt); 
   
   // Accessors.
   double C() const { return mC; }
   double Rmax() const { return mRmax; }
   double rhoMin() const { return mRhoMin; }

   //! Equivalence.
   virtual bool operator==(const Spheral::UpdatePolicyBase<Dim<3>, FieldSpace::Field<Dim<3>, Dim<3>::Scalar> >& rhs) const;

   private:

   // Parameters.
   double mC, mRmax, mRhoMin;

   // Disabled methods.
   SpitzerResistivityUpdatePolicy(const SpitzerResistivityUpdatePolicy& rhs);
   SpitzerResistivityUpdatePolicy& operator=(const SpitzerResistivityUpdatePolicy& rhs);
};

}
}

#else

// Forward declaration.
namespace Spheral {
namespace MHDSpace {
   class SpitzerResistivityUpdatePolicy;
}
}

#endif
