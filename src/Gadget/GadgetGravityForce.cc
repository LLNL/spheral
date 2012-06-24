//---------------------------------Spheral++----------------------------------//
// GravityForce -- Gadget's gravitational interactions.
// implementations.
//
//! \author $Author: mikeowen $
//! \version $Revision: 622 $
//! \date $Date: 2003-08-31 21:14:54 -0700 (Sun, 31 Aug 2003) $
//----------------------------------------------------------------------------//

#ifdef USE_GADGET

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include "GadgetGravityForce.hh"
#include "FakeGadgetICFile.hh"
#include "Geometry/Dimension.hh"
#include "DataBase.hh"
#include "FieldList.hh"
#include "Field.hh"
#include "cdebug.hh"
#include "DBC.hh"

extern "C" {
#include "allvars.h" // Gadget's global variables.
#include "proto.h"   // Gadget's function prototypes.
}

//------------------------------------------------------------------------------
//                          Implementation notes for Gadget
//------------------------------------------------------------------------------
// 1. Gadget uses cgs units, so Spheral needs to be in cgs mode in order 
//    to work properly.
// 2. Gadget should be built with the boundary conditions corresponding to 
//    the test problem to be run.  The types of BC's are free and periodic.
// 3. Do NOT compile Gadget with the -DDIAG option!  This causes a 
//    segmentation violation during the printing of diagnostics in the 
//    calculation of the gravity tree!
//------------------------------------------------------------------------------

namespace Spheral {
namespace GadgetSpace {

bool GravityForce::mIsInitialized = false;

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
GravityForce::GravityForce():
  mSizeOfGadgetParticleList(0),
  mPotential(FieldList<Dim<3>, Scalar>::Copy),
  mExtraEnergy(0.0) {
  

  using namespace std;
  cdebug << "GravityForce::GravityForce()" << endl;

  // NOTE: We will not initialize Gadget until the first time that we 
  // use it to evaluate our derivatives.  Originally, we had intended to 
  // initialize it with a zero-length particle list, but Gadget has a fit 
  // when you do this, so we were forced to defer initialization until 
  // its first use.  Oh, well.
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
GravityForce::
~GravityForce() {
  cdebug << "GravityForce::~GravityForce()" << endl;
  // Delete the existing particles from Gadget.
  for (int i = 0; i < NumPart; ++i)
  {
    if (P != 0)
    {
      delete [] P;
    } // end if
  } // end for

  // No more Gadget object!
  mIsInitialized = false;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void 
GravityForce::
evaluateDerivatives(const DataBase<Dim<3> >& dataBase,
                    const Scalar& time,
                    const Scalar& dt,
                    FieldList<Dim<3>, Vector>& DvDt) const
{
  cdebug << "Gadget: Evaluating derivatives...\n";
  // Important variables in Gadget:
  // P - array of particle data
  // NumPart - number of particles (per process)
  // Translate each particle in Spheral to a particle in Gadget.
  int i = 0;
  const FieldList<Dim<3>, Vector>& position = 
    dataBase.globalPosition();
  const FieldList<Dim<3>, Scalar>& mass = dataBase.globalMass();
  cdebug << "Gadget: Translating particle quantities...\n";
  for (DataBase<Dim<3> >::IDIterator 
       nodeIter = dataBase.internalNodeBegin();
       nodeIter != dataBase.internalNodeEnd();
       ++nodeIter) {
    // Position.
    CHECK(dataBase.globalPosition().numInternalNodes() == NumPart);
    for (int j = 0; j < Dim<3>::nDim; ++j) {
      // FIXME: Unit conversion here!
      P[i].Pos[j] = static_cast<float>(position(nodeIter)(j));
    } // end for

    // Mass.  FIXME: Does the mass of a spheral particle change under 
    // normal circumstances?
    P[i].Mass = static_cast<float>(mass(nodeIter));
    ++i;
  } // end for

  // Compute the gravitational forces of the particles upon one another.
  cdebug << "Gadget: Computing particle accelerations...\n";
  find_next_time();
  compute_accelerations(0);

  // Kerplunk the particle accelerations into Spheral's DvDt, and the 
  // particle gravitational potentials into mPotential.
  i = 0;
  cdebug << "Gadget: Imparting accelarations to Spheral particles...\n";
  for (DataBase<Dim<3> >::IDIterator 
       nodeIter = dataBase.nodeBegin();
       nodeIter != dataBase.nodeEnd();
       ++nodeIter) {

    // FIXME: Convert units back!
    
    // Accelerations.
    for (int j = 0; j < Dim<3>::nDim; ++j)
    {
      (DvDt(nodeIter))(j) = static_cast<double>(P[i].Accel[j]);
    } // end for

    // Gravitational potential.
    mPotential(nodeIter) = static_cast<double>(P[i].Potential);
    ++i;
  } // end for

  // Compute the potential for each of Gadget's particles.
  cdebug << "Computing potential...\n";
  compute_potential();

  // Compute the "global quantities" of Gadget, including the total 
  // potential energy.
  cdebug << "Gadget: Computing global quantities...\n";
  compute_global_quantities_of_system();
  // FIXME: Need units conversion here, as well!
  mExtraEnergy = static_cast<double>(SysState.EnergyPot);
  // We're done!
}

//------------------------------------------------------------------------------
// initialize()
//------------------------------------------------------------------------------
void 
GravityForce::
initialize(const DataBase<Dim<3> >& db, 
           ConstBoundaryIterator boundaryBegin,
           ConstBoundaryIterator boundaryEnd,
           const Scalar& time, const Scalar& dt)
{
  cdebug << "Gadget: Initializing...\n";
  // If we haven't yet initialized Gadget, we should do it now.
  if (!mIsInitialized) 
  {
    mInitializeGadget(db);
  } // end if
  else
  {
    // Ensure that Gadget's particle list is sized properly.  In particular, 
    // we can use the existing particle list if it is large enough to accomodate
    // all of Spheral's particles.  If Spheral has more particles than Gadget,
    // we have to resize Gadget's particle list so that it is large enough.
    // FIXME: Is there an easier way to get at the total number of nodes in 
    // Spheral???
    size_t numberOfInternalSpheralNodes = 0;
    for (DataBase<Dim<3> >::ConstNodeListIterator 
        nodeListIter = db.nodeListBegin();
        nodeListIter != db.nodeListEnd();
        ++nodeListIter) {
      numberOfInternalSpheralNodes += (*nodeListIter)->numInternalNodes();
    } // end for

    if (mSizeOfGadgetParticleList < numberOfInternalSpheralNodes)
    {
      if (P != 0)
      {
        free(P);
        P = 0;
      } // end if
      mSizeOfGadgetParticleList = numberOfInternalSpheralNodes;
      P = (particle_data*)malloc(sizeof(particle_data)*mSizeOfGadgetParticleList);

      // ForceFlag and CurrentTime have to be set properly for every particle.
      for (int i = 0; i < mSizeOfGadgetParticleList; ++i)
      {
        P[i].ForceFlag = i + 1;
        P[i].CurrentTime = All.TimeBegin;
      } // end for
    } // end if
    NumPart = static_cast<int>(numberOfInternalSpheralNodes);
  } // end if

  // Resize the gravitational potential if necessary.  This can be done 
  // by setting the fieldlist's nodelists.
  if (mPotential.numFields() == 0)
  {
    for (DataBase<Dim<3> >::ConstNodeListIterator 
         nodeListIter = db.nodeListBegin();
         nodeListIter != db.nodeListEnd();
         ++nodeListIter)
    {
      Spheral::FieldSpace::Field<Dim<3>, Scalar> potentialField(**nodeListIter);
      mPotential.appendField(potentialField);
    } // end for
  } // end if
  cdebug << "Gadget: Done initializing\n";
}

//------------------------------------------------------------------------------
// dt()
//------------------------------------------------------------------------------
GravityForce::Scalar 
GravityForce::
dt(const DataBase<Dim<3> >& dataBase, Scalar currentTime) const {
  // Let's try using the dynamical timestep, since Gadget appears not to put 
  // any computed upper bound on anything.
  // The dynamical time is dt ~ 1/sqrt(G*rhoMax), where rhoMax is the 
  // maximum mass density of the particle distribution.  We can nab this 
  // from Spheral, and just use G in cgs.

  // What the heck is the maximum mass density?  Let's find out.  We have to 
  // go through the somewhat arduous process of going through all the fields 
  // in the massDensity fieldlist, getting those maximum values, and then 
  // finding the maximum of the maximums.  Stupid FieldList nonstandard
  // containers.  Grumble grumble.
  FieldList<Dim<3>, Scalar> massDensity = dataBase.fluidMassDensity();
  double rhoMax = 0.0;
  for (FieldList<Dim<3>, Scalar>::iterator 
       fieldPtrIter = massDensity.begin();
       fieldPtrIter != massDensity.end();
       ++fieldPtrIter)
  {
    double maxForThisField = 
      *(std::max_element((*fieldPtrIter)->begin(),
                       (*fieldPtrIter)->end()));
    if (rhoMax < maxForThisField)
    {
      rhoMax = maxForThisField;
    } // end if
  } // end for
  double G = static_cast<double>(All.G);
  return 1.0/std::sqrt(G * rhoMax);
}

//------------------------------------------------------------------------------
// extraEnergy()
//------------------------------------------------------------------------------
GravityForce::Scalar 
GravityForce::
extraEnergy() const
{
  return mExtraEnergy;
}

//------------------------------------------------------------------------------
// extraEnergy()
//------------------------------------------------------------------------------
const FieldList<Dim<3>, GravityForce::Scalar>&
GravityForce::
potential() const
{
  return mPotential;
}

//------------------------------------------------------------------------------
// valid()
//------------------------------------------------------------------------------
bool 
GravityForce::
valid() const
{
  // We're always valid, sir!  (This is crap, but we can make no other 
  // assumptions right now.)
  return true;
}

//------------------------------------------------------------------------------
void
GravityForce::
mInitializeGadget(const DataBase<Dim<3> >& db)
{
  // Initialize Gadget.
  //
  // We have to fake an input file and feed it in 
  // to properly initialize all of Gadget's innards.  
  // The FakeGadgetICFile object creates a bogus IC file and then 
  // deletes this file upon its destruction.  
  FakeGadgetICFile fakeIC(db);

  // Now tweak Gadget's parameters and call it's init() function.
  //
  // Gadget should know what IC format to use.
  All.ICFormat = 1;

  // Gadget should know the name of the bogus file.
  strncpy(All.InitCondFile, fakeIC.name(), 100);
  
  // For serial Gadget, the particle allocation factor is 1. 
  All.PartAllocFactor = 1.0;

  // Error tolerances.  Cribbed from Gadget/.../parameterfiles/galaxy.tex.
  All.ErrTolIntAccuracy = 1.0;
  All.ErrTolVelScale = 10.0;
  All.MaxSizeTimestep = 0.01;
  All.ErrTolTheta = 0.8;
  All.TypeOfOpeningCriterion = 1;
  All.ErrTolForceAcc = 0.02;
  All.MaxNodeMove = 0.05;
  All.TreeUpdateFrequency = 0.1;

  // Gadget's method of timestep selection.
  All.TypeOfTimestepCriterion = 1;
  
  // The tree allocation factor should be about 0.7, according to the 
  // documentation.  However, our small-particle-distribution tests fail 
  // to allocate enough space for trees, so we've bumped it up to 4.0 for 
  // safety.
  All.TreeAllocFactor = 4.0;

  // Turn off comoving integration.
  All.ComovingIntegrationOn = 0;

  // Softening scale for Spheral particles.
  All.SofteningHalo = 1.0;
  All.SofteningHaloMaxPhys = 1.0;

  // Set units to cgs.
  // FIXME: If we want to give Gadget other units, we do it here!
  All.UnitMass_in_g = 1.0;
  All.UnitTime_in_s = 1.0;
  All.UnitLength_in_cm = 1.0;
  All.UnitVelocity_in_cm_per_s = 1.0;

  // Hubble parameters. Probably incorrect.
  All.Hubble = 1.0;
  All.HubbleParam = 1.0;
  All.Omega0 = 0.0;

  // Tell Gadget to read masses from our IC's.
  All.MassTable[0] = 0;

  // Call Gadget's init().  
  init();

  // We have initialized the Gadget object.
  mIsInitialized = true;
  mSizeOfGadgetParticleList = 0;
  for (DataBase<Dim<3> >::ConstNodeListIterator 
      nodeListIter = db.nodeListBegin();
      nodeListIter != db.nodeListEnd();
      ++nodeListIter) {
    mSizeOfGadgetParticleList += (*nodeListIter)->numInternalNodes();
  } // end for
  
  ENSURE(mIsInitialized);
}
//------------------------------------------------------------------------------

} // end namespace GadgetSpace
} // end namespace Spheral

//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------


#endif // USE_GADGET
