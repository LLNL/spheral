//---------------------------------Spheral++----------------------------------//
// MFVHydroBase -- spheralized verions of "Meshless Finite Mass" 
//   Hopkins P.F. (2015) "A New Class of Accurate, Mesh-Free Hydrodynamic 
//   Simulation Methods," MNRAS, 450(1):53-110
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"
#include "NodeList/SmoothingScaleBase.hh"
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

#include "GSPH/MFVHydroBase.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "GSPH/computeSumVolume.hh"
#include "GSPH/computeMFMDensity.hh"
#include "GSPH/Policies/ReplaceWithRatioPolicy.hh"
#include "GSPH/Policies/IncrementSpecificFromTotalPolicy.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <sstream>

using std::string;
using std::min;
using std::max;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MFVHydroBase<Dimension>::
MFVHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
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
  mDmassDt(FieldStorageType::CopyFields),
  mDthermalEnergyDt(FieldStorageType::CopyFields),
  mDmomentumDt(FieldStorageType::CopyFields),
  mDvolumeDt(FieldStorageType::CopyFields){
    mDmassDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::mass);
    mDthermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + GSPHFieldNames::thermalEnergy);
    mDmomentumDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + GSPHFieldNames::momentum);
    mDvolumeDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::volume);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
MFVHydroBase<Dimension>::
~MFVHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  GenericRiemannHydro<Dimension>::initializeProblemStartup(dataBase);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  GenericRiemannHydro<Dimension>::registerState(dataBase,state);
  
  auto massDensity = dataBase.fluidMassDensity();
  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);

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

  auto rhoPolicy = std::make_shared<ReplaceWithRatioPolicy<Dimension,Scalar>>(HydroFieldNames::mass,
                                                                         HydroFieldNames::volume,
                                                                         HydroFieldNames::volume);

  auto massPolicy = std::make_shared<IncrementFieldList<Dimension, Vector>>(true);
  auto momentumPolicy = std::make_shared<IncrementSpecificFromTotalPolicy<Dimension, Vector>>(HydroFieldNames::velocity, IncrementFieldList<Dimension, Vector>::prefix() + GSPHFieldNames::momentum);
  auto thermalEnergyPolicy = std::make_shared<IncrementSpecificFromTotalPolicy<Dimension, Scalar>>(HydroFieldNames::specificThermalEnergy, IncrementFieldList<Dimension, Scalar>::prefix() + GSPHFieldNames::thermalEnergy);

  // special policies to get velocity update from momentum deriv
  state.enroll(GSPHFieldNames::momentumPolicy,momentumPolicy);
  state.enroll(GSPHFieldNames::thermalEnergyPolicy,thermalEnergyPolicy);

  // normal state variables
  state.enroll(mass,        massPolicy);
  state.enroll(massDensity, rhoPolicy);
  state.enroll(volume,      volumePolicy);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::registerDerivatives(dataBase,derivs);
  dataBase.resizeFluidFieldList(mDmassDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::mass, false);
  dataBase.resizeFluidFieldList(mDthermalEnergyDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + GSPHFieldNames::thermalEnergy, false);
  dataBase.resizeFluidFieldList(mDmomentumDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + GSPHFieldNames::momentum, false);
  dataBase.resizeFluidFieldList(mDvolumeDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::volume, false);
  derivs.enroll(mDmassDt);
  derivs.enroll(mDthermalEnergyDt);
  derivs.enroll(mDmomentumDt);
  derivs.enroll(mDvolumeDt);
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
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
MFVHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::initialize(time,dt,dataBase,state,derivs); 
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
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
MFVHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::applyGhostBoundaries(state,derivs);
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  GenericRiemannHydro<Dimension>::enforceBoundaries(state,derivs);
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  GenericRiemannHydro<Dimension>::dumpState(file,pathName);
  file.write(mDmassDt, pathName + "/DmassDt");
  file.write(mDthermalEnergyDt, pathName + "/DthermalEnergyDt");
  file.write(mDmomentumDt, pathName + "/DmomentumDt");
  file.write(mDvolumeDt, pathName + "/DvolumeDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  GenericRiemannHydro<Dimension>::restoreState(file,pathName);
  file.read(mDmassDt, pathName + "/DmassDt");
  file.read(mDthermalEnergyDt, pathName + "/DthermalEnergyDt");
  file.read(mDmomentumDt, pathName + "/DmomentumDt");
  file.read(mDvolumeDt, pathName + "/DvolumeDt");
}

}
