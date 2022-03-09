//---------------------------------Spheral++----------------------------------//
// MFMHydroBase -- The SPH/ASPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "SPH/computeSPHSumMassDensity.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"

#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"

#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/Timer.hh"

#include "GSPH/MFMHydroBase.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "GSPH/computeSumVolume.hh"
#include "GSPH/computeMFMDensity.hh"
#include "GSPH/Policies/ReplaceWithRatioPolicy.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>

using std::vector;
using std::string;
using std::pair;
using std::to_string;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MFMHydroBase<Dimension>::
MFMHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
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
                                 densityUpdate,
                                 HUpdate,
                                 epsTensile,
                                 nTensile,
                                 xmin,
                                 xmax),
  mDvolumeDt(FieldStorageType::CopyFields){
    
    mDvolumeDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::volume);


}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
MFMHydroBase<Dimension>::
~MFMHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  GenericRiemannHydro<Dimension>::initializeProblemStartup(dataBase);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  GenericRiemannHydro<Dimension>::registerState(dataBase,state);
  
  auto massDensity = dataBase.fluidMassDensity();
  auto volume = state.fields(HydroFieldNames::volume, 0.0);

  std::shared_ptr<CompositeFieldListPolicy<Dimension, Scalar> > volumePolicy(new CompositeFieldListPolicy<Dimension, Scalar>());
  for (auto itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    auto massi = (*itr)->mass();
    auto minVolume = massi.min()/(*itr)->rhoMax();
    auto maxVolume = massi.max()/(*itr)->rhoMin();
    volumePolicy->push_back(new IncrementBoundedState<Dimension, Scalar>(minVolume,
                                                                         maxVolume));
  }

  PolicyPointer rhoPolicy(new ReplaceWithRatioPolicy<Dimension,Scalar>(HydroFieldNames::mass,
                                                                       HydroFieldNames::volume,
                                                                       HydroFieldNames::volume));
  
  // normal state variables
  state.enroll(massDensity, rhoPolicy);
  state.enroll(volume, volumePolicy);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::registerDerivatives(dataBase,derivs);

  dataBase.resizeFluidFieldList(mDvolumeDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::volume, false);
  derivs.enroll(mDvolumeDt);
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::preStepInitialize(dataBase,state,derivs);

  if(this->densityUpdate() == MassDensityType::RigorousSumDensity){
    // plop into an intialize volume function
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
          auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
          auto  volume = state.fields(HydroFieldNames::volume, 0.0);

    computeSumVolume(dataBase.connectivityMap(),this->kernel(),position,H,volume);
    computeMFMDensity(mass,volume,massDensity);
  
    for (auto boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr){
      (*boundaryItr)->applyFieldListGhostBoundary(volume);
      (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    }
    for (auto boundaryItr = this->boundaryBegin(); 
         boundaryItr < this->boundaryEnd(); 
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Initialize the hydro before calling evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::initialize(time,dt,dataBase,state,derivs); 
    // plop into an intialize volume function
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
          auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
          auto  volume = state.fields(HydroFieldNames::volume, 0.0);

    computeSumVolume(dataBase.connectivityMap(),this->kernel(),position,H,volume);
    computeMFMDensity(mass,volume,massDensity);
  
    for (auto boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr){
      (*boundaryItr)->applyFieldListGhostBoundary(volume);
      (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    }
    for (auto boundaryItr = this->boundaryBegin(); 
         boundaryItr < this->boundaryEnd(); 
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  GenericRiemannHydro<Dimension>::finalizeDerivatives(time,dt,dataBase,state,derivs);
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::applyGhostBoundaries(state,derivs);
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::enforceBoundaries(state,derivs);
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  GenericRiemannHydro<Dimension>::dumpState(file,pathName);

  file.write(mDvolumeDt, pathName + "/DvolumeDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFMHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  GenericRiemannHydro<Dimension>::restoreState(file,pathName);

  file.read(mDvolumeDt, pathName + "/DvolumeDt");
}

}

