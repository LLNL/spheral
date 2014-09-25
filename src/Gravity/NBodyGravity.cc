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
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "Utilities/DBC.hh"
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
NBodyGravity(const double plummerSofteningLength,
             const double maxDeltaVelocity,
             const double G):
  mPotential(FieldSpace::Copy),
  mExtraEnergy(0.0),
  mMaxDeltaVelocityFactor(maxDeltaVelocity),
  mSofteningLength(plummerSofteningLength),
  mG(G) {
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
                    StateDerivatives<Dimension>& derivs) const {
  using namespace NodeSpace;

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
  const unsigned numNodeLists = dataBase.numNodeLists();
  unsigned ifield = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator iitr = dataBase.fluidNodeListBegin();
       iitr != dataBase.fluidNodeListEnd();
       ++iitr, ++ifield) {
    const NodeList<Dimension>& nodeListi = **iitr;
    const unsigned firstGhostNodei = nodeListi.firstGhostNode();

    for (unsigned i = 0; i != firstGhostNodei; ++i) {
      unsigned jfield = 0;

      // Set the position derivative.
      DxDt(ifield, i) += velocity(ifield, i);

      // Get a reference to the acceleration vector.
      Vector& acceleration = DvDt(ifield, i);
    
      // Zero out the potential of this particle.
      mPotential(ifield, i) = 0.0;

      for (typename DataBase<Dimension>::ConstFluidNodeListIterator jitr = dataBase.fluidNodeListBegin();
           jitr != dataBase.fluidNodeListEnd();
           ++jitr, ++jfield) {
        const NodeList<Dimension>& nodeListj = **jitr;
        const unsigned nj = nodeListj.numNodes();
        for (unsigned j = 0; j != nj; ++j) {

          // Particles can't self-interact, silly!
          if (ifield != jfield or i != j) {

            // Contribute to the acceleration of this particle.
            const Vector r = position(ifield, i) - position(jfield, j);

            CHECK(r.magnitude2() != 0.0);
            Vector rHat = r.unitVector();
            Scalar distance2 = r.magnitude2() + softeningLength2;
            CHECK(distance2 != 0.0);
            acceleration -= mG * mass(jfield, j) * rHat / distance2;

            // Also sum up contributions to the potential and 
            // total potential energy.
            mPotential(ifield, i) -= mG * mass(jfield, j) / std::sqrt(distance2);
          }
        }
      }

      mExtraEnergy += mG * mass(ifield, i) * mPotential(ifield, i);

      // Capture the maximum acceleration and velocity magnitudes.
      const Scalar accelMagnitude = acceleration.magnitude();
      if (mOldMaxAcceleration < accelMagnitude) {
        mOldMaxAcceleration = accelMagnitude + 1e-10;
        CHECK(mOldMaxAcceleration != 0.0);
      } // end if

      const Scalar velocityMagnitude = velocity(ifield, i).magnitude();
      if (mOldMaxVelocity < velocityMagnitude) {
        mOldMaxVelocity = velocityMagnitude;
      }
    }
  }
}

//------------------------------------------------------------------------------
// initialize()
//------------------------------------------------------------------------------
template <typename Dimension>
void 
NBodyGravity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& db) {

  // Allocate space for the gravitational potential FieldList.
  mPotential = db.newGlobalFieldList(0.0, "gravitational potential");

}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
typename NBodyGravity<Dimension>::TimeStepType
NBodyGravity<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {

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
extraEnergy() const {
  return mExtraEnergy;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
const FieldList<Dimension, typename NBodyGravity<Dimension>::Scalar>&
NBodyGravity<Dimension>::
potential() const {
  return mPotential;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
bool 
NBodyGravity<Dimension>::
valid() const {
  // We're always valid, sir!  (This is crap, but we can make no other 
  // assumptions right now.)
  return true;
}

//------------------------------------------------------------------------------
template <typename Dimension>
double
NBodyGravity<Dimension>::
G() const {
  return mG;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
double
NBodyGravity<Dimension>::
softeningLength() const {
  return mSofteningLength;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
void
NBodyGravity<Dimension>::
softeningLength(const double x) {
  VERIFY(x >= 0.0);
  mSofteningLength = x;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
} // end namespace GravitySpace
} // end namespace Spheral

