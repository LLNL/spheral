//---------------------------------Spheral++----------------------------------//
// NBodyGravity implementation.
//
//! \author $Author: mikeowen $
//! \version $Revision: 4012 $
//! \date $Date: 2011-02-20 10:58:32 -0800 (Sun, 20 Feb 2011) $
//----------------------------------------------------------------------------//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "NBodyGravity.hh"
#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "cdebug.hh"
#include "DBC.hh"
#include "Material/PhysicalConstants.hh"

namespace Spheral {
namespace GravitySpace {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
template <typename Dimension>
NBodyGravity<Dimension>::
NBodyGravity(double plummerSofteningLength,
             double maxDeltaVelocity):
  mPotential(FieldList<Dimension, Scalar>::Copy),
  mExtraEnergy(0.0),
  mMaxDeltaVelocityFactor(maxDeltaVelocity),
  mSofteningLength(plummerSofteningLength) {
  using namespace std;
  cdebug << "NBodyGravity<Dimension>::NBodyGravity()" << endl;

}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
NBodyGravity<Dimension>::
~NBodyGravity() {
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
void 
NBodyGravity<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const
{
  using namespace NodeSpace;

  // Don't forget Cavendish's constant.
  // FIXME: Currently in cgs units.
  Scalar G = 6.672e-8;

  // Find the square of the Plummer softening length.
  Scalar softeningLength2 = mSofteningLength * mSofteningLength;

  // Access to pertinent fields in the database.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);

  // Get the acceleration and position change vectors we'll be modifying.
  FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Zero out the total gravitational potential energy.
  mExtraEnergy = 0.0;

  // Loop over each particle...
  for (InternalNodeIterator<Dimension> ithNodeIter = dataBase.internalNodeBegin();
       ithNodeIter != dataBase.internalNodeEnd();
       ++ithNodeIter)
  {

    // Set the position derivative.
    DxDt(ithNodeIter) += velocity(ithNodeIter);

    // Get a reference to the acceleration vector.
    Vector& acceleration = DvDt(ithNodeIter);
    
    // Zero out the potential of this particle.
    mPotential(ithNodeIter) = 0.0;

    // Now loop over every particle that's not this one...
    for (InternalNodeIterator<Dimension> jthNodeIter = dataBase.internalNodeBegin();
         jthNodeIter != dataBase.internalNodeEnd();
         ++jthNodeIter)
    {
      // Particles can't self-interact, silly!
      if (ithNodeIter != jthNodeIter)
      {
        // Contribute to the acceleration of this particle.
        Vector r = position(ithNodeIter) - position(jthNodeIter);

        CHECK(r.magnitude2() != 0.0);
        Vector rHat = r.unitVector();
        Scalar distance2 = r.magnitude2() + softeningLength2;
        CHECK(distance2 != 0.0);
        acceleration -= G * mass(jthNodeIter) * rHat / distance2;

        // Also sum up contributions to the potential and 
        // total potential energy.
        mPotential(ithNodeIter) += G * mass(jthNodeIter) / std::sqrt(distance2);
        mExtraEnergy += G * mass(ithNodeIter) * mPotential(ithNodeIter);
      } // end if
    } // end for

    // Capture the maximum acceleration and velocity magnitudes.
    Scalar accelMagnitude = acceleration.magnitude();
    if (mOldMaxAcceleration < accelMagnitude)
    {
      mOldMaxAcceleration = accelMagnitude + 1e-10;
      CHECK(mOldMaxAcceleration != 0.0);
    } // end if
    Scalar velocityMagnitude = velocity(ithNodeIter).magnitude();
    if (mOldMaxVelocity < velocityMagnitude)
    {
      mOldMaxVelocity = velocityMagnitude;
    } // end if
  } // end for
}

//------------------------------------------------------------------------------
// initialize()
//------------------------------------------------------------------------------
template <typename Dimension>
void 
NBodyGravity<Dimension>::
initialize(const Scalar& time, 
           const Scalar& dt,
           const DataBase<Dimension>& db, 
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs)
{
  // Allocate space for the gravitational potential FieldList if necessary.
  if (mPotential.numFields() == 0)
  {
    mPotential = db.newGlobalFieldList(Scalar());
  } // end if
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
typename NBodyGravity<Dimension>::TimeStepType
NBodyGravity<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const
{

  // The maximum change in our velocity for the next time cycle is given 
  // by the mMaxDeltaVelocityFactor plus the max velocity and the max 
  // accelation from the last cycle.
  const double deltat = (mOldMaxVelocity*mOldMaxAcceleration/(mOldMaxAcceleration*mOldMaxAcceleration + 1.0e-10)) * mMaxDeltaVelocityFactor;

  stringstream reasonStream;
  reasonStream << "velocity: " << mOldMaxVelocity
               << ", acceleration: " << mOldMaxAcceleration
               << "dt = f*v/a: " << deltat
               << ends;
  return TimeStepType(deltat, reasonStream.str());
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
typename NBodyGravity<Dimension>::Scalar 
NBodyGravity<Dimension>::
extraEnergy() const
{
  return mExtraEnergy;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
const FieldList<Dimension, typename NBodyGravity<Dimension>::Scalar>&
NBodyGravity<Dimension>::
potential() const
{
  return mPotential;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
bool 
NBodyGravity<Dimension>::
valid() const
{
  // We're always valid, sir!  (This is crap, but we can make no other 
  // assumptions right now.)
  return true;
}

//------------------------------------------------------------------------------
} // end namespace GravitySpace
} // end namespace Spheral

//------------------------------------------------------------------------------
// Explict instantiations.
//------------------------------------------------------------------------------

namespace Spheral {
  namespace GravitySpace {

    template class NBodyGravity<Dim<1> >;
    template class NBodyGravity<Dim<2> >;
    template class NBodyGravity<Dim<3> >;

  } // end namespace GravitySpace
} // end namespace Spheral

