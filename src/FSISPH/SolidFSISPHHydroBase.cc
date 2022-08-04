//---------------------------------Spheral++----------------------------------//
// SolidFSISPHHydroBase -- SolidSPHHydro modified for large density discontinuities
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"

#include "SPH/SPHHydroBase.hh"                 
#include "SPH/SolidSPHHydroBase.hh"

#include "GSPH/Policies/PureReplaceFieldList.hh"
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "SolidMaterial/SolidEquationOfState.hh" 

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"

#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/timerLayer.hh"

#include "FSISPH/SolidFSISPHHydroBase.hh"
#include "FSISPH/FSIFieldNames.hh"
#include "FSISPH/computeFSISPHSumMassDensity.hh"
#include "FSISPH/computeHWeightedFSISPHSumMassDensity.hh"
#include "FSISPH/SlideSurface.hh"


#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;


namespace Spheral {


inline
Dim<1>::SymTensor
tensileStressCorrection(const Dim<1>::SymTensor& sigma) {
  if (sigma.xx() > 0.0) {
    return -sigma;
  } else {
    return Dim<1>::SymTensor::zero;
  }
}

inline
Dim<2>::SymTensor
tensileStressCorrection(const Dim<2>::SymTensor& sigma) {
  const EigenStruct<2> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  Dim<2>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,
                           0.0,                              (lambday > 0.0 ? -lambday : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

inline
Dim<3>::SymTensor
tensileStressCorrection(const Dim<3>::SymTensor& sigma) {
  const EigenStruct<3> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  const double lambdaz = eigen.eigenValues.z();
  Dim<3>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,                              0.0,
                           0.0,                              (lambday > 0.0 ? -lambday : 0.0), 0.0,
                           0.0,                              0.0,                              (lambdaz > 0.0 ? -lambdaz : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}


//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SolidFSISPHHydroBase<Dimension>::
SolidFSISPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                  DataBase<Dimension>& dataBase,
                  ArtificialViscosity<Dimension>& Q,
                  SlideSurface<Dimension>& slides,
                  const TableKernel<Dimension>& W,
                  const double filter,
                  const double cfl,
                  const double surfaceForceCoefficient,
                  const double densityStabilizationCoefficient,
                  const double specificThermalEnergyDiffusionCoefficient,
                  const double xsphCoefficient,
                  const InterfaceMethod interfaceMethod,
                  const KernelAveragingMethod kernelAveragingMethod,
                  const std::vector<int> sumDensityNodeLists,
                  const bool useVelocityMagnitudeForDt,
                  const bool compatibleEnergyEvolution,
                  const bool evolveTotalEnergy,
                  const bool gradhCorrection,
                  const bool XSPH,
                  const bool correctVelocityGradient,
                  const MassDensityType densityUpdate,
                  const HEvolutionType HUpdate,
                  const double epsTensile,
                  const double nTensile,
                  const bool damageRelieveRubble,
                  const bool strengthInDamage,
                  const Vector& xmin,
                  const Vector& xmax):
  SolidSPHHydroBase<Dimension>(smoothingScaleMethod,
                               dataBase,
                               Q,
                               W,
                               W,  //WPi
                               W,  //WGrad
                               filter,
                               cfl,
                               useVelocityMagnitudeForDt,
                               compatibleEnergyEvolution,
                               evolveTotalEnergy,
                               gradhCorrection, 
                               XSPH,
                               correctVelocityGradient,
                               true, // sumMassDensityOverAllNodeLists
                               densityUpdate,
                               HUpdate,
                               epsTensile,
                               nTensile,
                               damageRelieveRubble,
                               strengthInDamage,
                               xmin,
                               xmax),
  mSlideSurface(slides),
  mSurfaceForceCoefficient(surfaceForceCoefficient),
  mDensityStabilizationCoefficient(densityStabilizationCoefficient),
  mSpecificThermalEnergyDiffusionCoefficient(specificThermalEnergyDiffusionCoefficient),
  mXSPHCoefficient(xsphCoefficient),
  mInterfaceMethod(interfaceMethod),
  mKernelAveragingMethod(kernelAveragingMethod),
  mApplySelectDensitySum(false),
  mSumDensityNodeLists(sumDensityNodeLists),
  mPairDepsDt(),
  mRawPressure(FieldStorageType::CopyFields),
  mDPDx(FieldStorageType::CopyFields),
  mDepsDx(FieldStorageType::CopyFields),
  mInterfaceNormals(FieldStorageType::CopyFields),
  mInterfaceFraction(FieldStorageType::CopyFields),
  mInterfaceSmoothness(FieldStorageType::CopyFields),
  mNewInterfaceNormals(FieldStorageType::CopyFields),
  mSmoothedInterfaceNormals(FieldStorageType::CopyFields),
  mNewInterfaceFraction(FieldStorageType::CopyFields),
  mNewInterfaceSmoothness(FieldStorageType::CopyFields){

    mPairDepsDt.clear();

    // see if we're summing density for any nodelist
    auto numNodeLists = dataBase.numNodeLists();
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      if (sumDensityNodeLists[nodeListi]==1){
        mApplySelectDensitySum = true;
      } 
    }
    
    mRawPressure = dataBase.newFluidFieldList(0.0, FSIFieldNames::rawPressure);
    mDPDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::pressureGradient);
    mDepsDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::specificThermalEnergyGradient);
    mInterfaceNormals = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::interfaceNormals);
    mInterfaceFraction = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceFraction);
    mInterfaceSmoothness = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceSmoothness);
    mNewInterfaceNormals = dataBase.newFluidFieldList(Vector::one, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals);
    mSmoothedInterfaceNormals = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::smoothedInterfaceNormals);
    mNewInterfaceFraction = dataBase.newFluidFieldList(0.0, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceFraction);
    mNewInterfaceSmoothness = dataBase.newFluidFieldList(0.0, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness);
    
  }

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SolidFSISPHHydroBase<Dimension>::
~SolidFSISPHHydroBase() {
}

//------------------------------------------------------------------------------
// initialization on problem start up
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase){

  SolidSPHHydroBase<Dimension>::initializeProblemStartup(dataBase);

  mRawPressure+=this->pressure();
  
}

//------------------------------------------------------------------------------
// Register states
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("SolidFSISPHregisterState");

  SolidSPHHydroBase<Dimension>::registerState(dataBase,state);
  

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;
  
  dataBase.resizeFluidFieldList(mRawPressure, 0.0, FSIFieldNames::rawPressure, false);
  dataBase.resizeFluidFieldList(mInterfaceNormals, Vector::zero, FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mInterfaceFraction, 0.0, FSIFieldNames::interfaceFraction,false); 
  dataBase.resizeFluidFieldList(mInterfaceSmoothness, 0.0, FSIFieldNames::interfaceSmoothness,false);
   
  PolicyPointer rawPressurePolicy(new PressurePolicy<Dimension>());
  PolicyPointer interfaceNormalsPolicy(new PureReplaceFieldList<Dimension,Vector>(ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals));
  PolicyPointer interfaceFractionPolicy(new PureReplaceFieldList<Dimension,Scalar>(ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceFraction));
  PolicyPointer interfaceSmoothnessPolicy(new PureReplaceFieldList<Dimension,Scalar>(ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness));
  
  // Override the specific thermal energy policy if compatible
  if(this->compatibleEnergyEvolution()){
    auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
    CHECK(specificThermalEnergy.numFields() == dataBase.numFluidNodeLists());
    PolicyPointer epsPolicy(new CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>(dataBase));
    state.enroll(specificThermalEnergy, epsPolicy);
  }

  state.enroll(mRawPressure,rawPressurePolicy);
  state.enroll(mInterfaceNormals,interfaceNormalsPolicy); 
  state.enroll(mInterfaceFraction,interfaceFractionPolicy);
  state.enroll(mInterfaceSmoothness,interfaceSmoothnessPolicy); 

  TIME_END("SolidFSISPHregisterState");
}

//------------------------------------------------------------------------------
// Register Derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>&  dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidFSISPHregisterDerivs");

  // Call the ancestor method.
  SolidSPHHydroBase<Dimension>::registerDerivatives(dataBase, derivs);  
  
  // make sure we're tracking the right number of node lists
  dataBase.resizeFluidFieldList(mDPDx, Vector::zero, FSIFieldNames::pressureGradient, false);
  dataBase.resizeFluidFieldList(mDepsDx, Vector::zero, FSIFieldNames::specificThermalEnergyGradient, false);
  dataBase.resizeFluidFieldList(mNewInterfaceNormals, Vector::zero,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mSmoothedInterfaceNormals, Vector::zero,  FSIFieldNames::smoothedInterfaceNormals,false);
  dataBase.resizeFluidFieldList(mNewInterfaceFraction, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() +  FSIFieldNames::interfaceFraction,false); 
  dataBase.resizeFluidFieldList(mNewInterfaceSmoothness, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness,false);

  // enroll 
  derivs.enrollAny(HydroFieldNames::pairWork,  mPairDepsDt);
  derivs.enroll(mDPDx);
  derivs.enroll(mDepsDx);
  derivs.enroll(mNewInterfaceNormals);
  derivs.enroll(mSmoothedInterfaceNormals);
  derivs.enroll(mNewInterfaceFraction);
  derivs.enroll(mNewInterfaceSmoothness);

  TIME_END("SolidFSISPHregisterDerivs");
}

//------------------------------------------------------------------------------
// FSI specialized density summmation
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_BEGIN("SolidFSISPHpreStepInitialize");
  if (mApplySelectDensitySum){
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
      const auto& W = this->kernel();
            auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeHWeightedFSISPHSumMassDensity(connectivityMap, W, mSumDensityNodeLists, position, mass, H, massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
  TIME_END("SolidFSISPHpreStepInitialize");
}

//------------------------------------------------------------------------------
// FSI specialized of the initialize method
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidFSISPHinitialize");

  const TableKernel<Dimension>& W = this->kernel();

  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  Q.initialize(dataBase, 
               state,
               derivs,
               this->boundaryBegin(),
               this->boundaryEnd(),
               time, 
               dt,
               W);
  
  // We depend on the caller knowing to finalize the ghost boundaries!
  TIME_END("SolidFSISPHinitialize");
}

//------------------------------------------------------------------------------
// For compatible energy we need to apply the bc conditions to acceleration
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
finalizeDerivatives(const Scalar /*time*/, 
                    const Scalar  /*dt*/,
                    const DataBase<Dimension>&  /*dataBase*/, 
                    const State<Dimension>& /*state*/,
                          StateDerivatives<Dimension>&  derivs) const {
                            
  if (this->compatibleEnergyEvolution()) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) {
           (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    }
    
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }

} // finalize


//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  SolidSPHHydroBase<Dimension>::applyGhostBoundaries(state,derivs);

  FieldList<Dimension, Scalar> interfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);
  FieldList<Dimension, Scalar> rawPressure = state.fields(FSIFieldNames::rawPressure, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceFraction);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceNormals);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceSmoothness);
    (*boundaryItr)->applyFieldListGhostBoundary(rawPressure);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  SolidSPHHydroBase<Dimension>::enforceBoundaries(state,derivs);

  FieldList<Dimension, Scalar> interfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);
  FieldList<Dimension, Scalar> rawPressure = state.fields(FSIFieldNames::rawPressure, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(interfaceFraction);
    (*boundaryItr)->enforceFieldListBoundary(interfaceNormals);
    (*boundaryItr)->enforceFieldListBoundary(interfaceSmoothness);
    (*boundaryItr)->enforceFieldListBoundary(rawPressure);
  }

}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  SolidSPHHydroBase<Dimension>::dumpState(file, pathName);

  file.write(mRawPressure, pathName + "/rawEosPressure");
  file.write(mDPDx, pathName + "/DpDx");
  file.write(mDepsDx, pathName + "/DepsDx");
  file.write(mInterfaceNormals, pathName + "/interfaceNormals");
  file.write(mInterfaceFraction, pathName + "/interfaceFraction");
  file.write(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  file.write(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.write(mSmoothedInterfaceNormals, pathName + "/smoothedInterfaceNormals");
  file.write(mNewInterfaceFraction, pathName + "/newInterfaceFraction");
  file.write(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  SolidSPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mRawPressure, pathName + "/rawEosPressure");
  file.read(mDPDx, pathName + "/DpDx");
  file.read(mDepsDx, pathName + "/DepsDx");
  file.read(mInterfaceNormals, pathName + "/interfaceNormals");
  file.read(mInterfaceFraction, pathName + "/interfaceFraction");
  file.read(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  file.read(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.read(mSmoothedInterfaceNormals, pathName + "/smoothedInterfaceNormals");
  file.read(mNewInterfaceFraction, pathName + "/newInterfaceFraction");
  file.read(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");
}

//------------------------------------------------------------------------------
// method for limited linear reconstruction between nodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
linearReconstruction(const typename Dimension::Vector& ri,
                     const typename Dimension::Vector& rj,
                     const typename Dimension::Scalar& yi,
                     const typename Dimension::Scalar& yj,
                     const typename Dimension::Vector& DyDxi,
                     const typename Dimension::Vector& DyDxj,
                           typename Dimension::Scalar& ytildei,
                           typename Dimension::Scalar& ytildej) const {
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto rij = (ri-rj);

  // relavant deltas in field value
  const auto Dy0 = (yi-yj);
  const auto Dyi = 0.5*DyDxi.dot(rij);
  const auto Dyj = 0.5*DyDxj.dot(rij);

  // ratios of SPH derivs to ij particle difference
  const auto denom = 2.0 / (sgn(Dy0) * std::max(tiny,abs(Dy0)));
  const auto xi = Dyi * denom;
  const auto xj = Dyj * denom;

  // limiter function - vanleer 1979
  const auto phii = ( xi > 0.0 ?  min(4.0*xi/((1.0 + xi)*(1.0 + xi)),1.0) : 0.0 );
  const auto phij = ( xj > 0.0 ?  min(4.0*xj/((1.0 + xj)*(1.0 + xj)),1.0) : 0.0 );        
  const auto phi = 0.5*(phii+phij);

  // linear constructed inteface values
  ytildei = yi - phi * Dyi;
  ytildej = yj + phi * Dyj;
}

} // Spheral namespace



