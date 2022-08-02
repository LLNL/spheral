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
#include "caliper/cali.h"

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
  CALI_MARK_BEGIN("GSPHinitializeStartup");
  GenericRiemannHydro<Dimension>::initializeProblemStartup(dataBase);
  CALI_MARK_END("GSPHinitializeStartup");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  CALI_MARK_BEGIN("GSPHregister");

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

  CALI_MARK_END("GSPHregister");
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  CALI_MARK_BEGIN("GSPHregisterDerivs");

  GenericRiemannHydro<Dimension>::registerDerivatives(dataBase,derivs);

  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  derivs.enroll(mDmassDensityDt);

  CALI_MARK_END("GSPHregisterDerivs");
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
  CALI_MARK_BEGIN("GSPHpreStepInitialize");
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
  CALI_MARK_END("GSPHpreStepInitialize");
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
  CALI_MARK_BEGIN("GSPHinitialize");

  GenericRiemannHydro<Dimension>::initialize(time,dt,dataBase,state,derivs);

  CALI_MARK_END("GSPHinitialize");
  
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
  CALI_MARK_BEGIN("GSPHfinalizeDerivs");
  GenericRiemannHydro<Dimension>::finalizeDerivatives(time,dt,dataBase,state,derivs);
  CALI_MARK_END("GSPHfinalizeDerivs");
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  CALI_MARK_BEGIN("GSPHghostBounds");

  GenericRiemannHydro<Dimension>::applyGhostBoundaries(state,derivs);

  CALI_MARK_END("GSPHghostBounds");
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  CALI_MARK_BEGIN("GSPHenforceBounds");

  GenericRiemannHydro<Dimension>::enforceBoundaries(state,derivs);

  CALI_MARK_END("GSPHenforceBounds");
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
