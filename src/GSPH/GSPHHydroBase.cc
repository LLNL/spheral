//---------------------------------Spheral++----------------------------------//
// GSPHHydroBase -- The Godunov SPH hydrodynamic package for Spheral++.
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "Hydro/HydroFieldNames.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/Timer.hh"

#include "GSPH/GSPHHydroBase.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "GSPH/computeSPHVolume.hh"
#include "GSPH/Policies/ReplaceWithRatioPolicy.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <sstream>

using std::string;
using std::min;
using std::max;

// Declare timers
extern Timer TIME_GSPH;
extern Timer TIME_GSPHinitializeStartup;
extern Timer TIME_GSPHregister;
extern Timer TIME_GSPHregisterDerivs;
extern Timer TIME_GSPHpreStepInitialize;
extern Timer TIME_GSPHinitialize;
extern Timer TIME_GSPHfinalizeDerivs;
extern Timer TIME_GSPHghostBounds;
extern Timer TIME_GSPHenforceBounds;
extern Timer TIME_GSPHevalDerivs;


namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GSPHHydroBase<Dimension>::
GSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
             DataBase<Dimension>& dataBase,
             RiemannSolverBase<Dimension>& riemannSolver,
             const TableKernel<Dimension>& W,
             const Scalar epsDiffusionCoeff,
             const double cfl,
             const bool useVelocityMagnitudeForDt,
             const bool compatibleEnergyEvolution,
             const bool evolveTotalEnergy,
             const bool XSPH,
             const bool correctVelocityGradient,
             const GradientType gradType,
             const MassDensityType densityUpdate,
             const HEvolutionType HUpdate,
             const double epsTensile,
             const double nTensile,
             const Vector& xmin,
             const Vector& xmax):
  GenericRiemannHydro<Dimension>(smoothingScaleMethod,
                                 dataBase,
                                 riemannSolver,
                                 W,
                                 epsDiffusionCoeff,
                                 cfl,
                                 useVelocityMagnitudeForDt,
                                 compatibleEnergyEvolution,
                                 evolveTotalEnergy,
                                 XSPH,
                                 correctVelocityGradient,
                                 gradType,
                                 densityUpdate,
                                 HUpdate,
                                 epsTensile,
                                 nTensile,
                                 xmin,
                                 xmax),
  mDmassDensityDt(FieldStorageType::CopyFields){

    mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GSPHHydroBase<Dimension>::
~GSPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  TIME_GSPHinitializeStartup.start();
  GenericRiemannHydro<Dimension>::initializeProblemStartup(dataBase);
  TIME_GSPHinitializeStartup.stop();
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_GSPHregister.start();

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  GenericRiemannHydro<Dimension>::registerState(dataBase,state);

  auto massDensity = dataBase.fluidMassDensity();
  auto volume = state.fields(HydroFieldNames::volume, 0.0);

  std::shared_ptr<CompositeFieldListPolicy<Dimension, Scalar> > rhoPolicy(new CompositeFieldListPolicy<Dimension, Scalar>());
  for (auto itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    rhoPolicy->push_back(new IncrementBoundedState<Dimension, Scalar>((*itr)->rhoMin(),
                                                                      (*itr)->rhoMax()));
  }

  PolicyPointer volumePolicy(new ReplaceWithRatioPolicy<Dimension,Scalar>(HydroFieldNames::mass,
                                                                        HydroFieldNames::massDensity,
                                                                        HydroFieldNames::massDensity));
  
  // normal state variables
  state.enroll(massDensity, rhoPolicy);
  state.enroll(volume, volumePolicy);

  TIME_GSPHregister.stop();
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_GSPHregisterDerivs.start();

  GenericRiemannHydro<Dimension>::registerDerivatives(dataBase,derivs);

  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  derivs.enroll(mDmassDensityDt);

  TIME_GSPHregisterDerivs.stop();
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  TIME_GSPHpreStepInitialize.start();
  GenericRiemannHydro<Dimension>::preStepInitialize(dataBase,state,derivs);

  if(this->densityUpdate() == MassDensityType::RigorousSumDensity){
    const auto& connectivityMap = dataBase.connectivityMap();
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
          auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
          auto  volume = state.fields(HydroFieldNames::volume, 0.0);

    computeSPHSumMassDensity(connectivityMap, this->kernel(), true, position, mass, H, massDensity);
    computeSPHVolume(mass,massDensity,volume);
  
    for (auto boundaryItr = this->boundaryBegin(); 
              boundaryItr < this->boundaryEnd(); 
              ++boundaryItr){
      (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      (*boundaryItr)->applyFieldListGhostBoundary(volume);
    }
    for (auto boundaryItr = this->boundaryBegin(); 
              boundaryItr < this->boundaryEnd(); 
              ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
  TIME_GSPHpreStepInitialize.stop();
}

//------------------------------------------------------------------------------
// Initialize the hydro before calling evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  TIME_GSPHinitialize.start();

  GenericRiemannHydro<Dimension>::initialize(time,dt,dataBase,state,derivs);

  TIME_GSPHinitialize.stop();
  
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_GSPHfinalizeDerivs.start();
  GenericRiemannHydro<Dimension>::finalizeDerivatives(time,dt,dataBase,state,derivs);
  TIME_GSPHfinalizeDerivs.stop();
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  TIME_GSPHghostBounds.start();

  GenericRiemannHydro<Dimension>::applyGhostBoundaries(state,derivs);

  TIME_GSPHghostBounds.stop();
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  TIME_GSPHenforceBounds.start();

  GenericRiemannHydro<Dimension>::enforceBoundaries(state,derivs);

  TIME_GSPHenforceBounds.stop();
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  GenericRiemannHydro<Dimension>::dumpState(file,pathName);

  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  GenericRiemannHydro<Dimension>::restoreState(file,pathName);

  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
}

}
