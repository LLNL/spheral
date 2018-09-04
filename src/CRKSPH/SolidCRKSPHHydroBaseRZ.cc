//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The CRKSPH/ACRKSPH solid material hydrodynamic
// package for Spheral++.
//
// This is the area-weighted RZ specialization.
//
// Created by JMO, Fri May 13 10:50:36 PDT 2016
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "CRKSPHHydroBase.hh"
#include "CRKSPHUtilities.hh"
#include "volumeSpacing.hh"
#include "computeHullVolumes.hh"
#include "computeCRKSPHSumVolume.hh"
#include "computeCRKSPHMoments.hh"
#include "computeCRKSPHCorrections.hh"
#include "computeCRKSPHSumMassDensity.hh"
#include "computeSolidCRKSPHSumMassDensity.hh"
#include "gradientCRKSPH.hh"
#include "Physics/GenericHydro.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/RZDeviatoricStressPolicy.hh"
#include "Strength/RZPlasticStrainPolicy.hh"
#include "Strength/BulkModulusPolicy.hh"
#include "Strength/ShearModulusPolicy.hh"
#include "Strength/YieldStrengthPolicy.hh"
#include "Strength/StrengthSoundSpeedPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "SPH/DamagedNodeCouplingWithFrags.hh"
#include "SolidMaterial/SolidEquationOfState.hh"

#include "SolidCRKSPHHydroBaseRZ.hh"

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

//------------------------------------------------------------------------------
// Compute the artificial tensile stress correction tensor for the given 
// stress tensor
//------------------------------------------------------------------------------
inline
Dim<1>::SymTensor
tensileStressCorrection(const Dim<1>::SymTensor& sigma) {
  if (sigma.xx() > 0.0) {
    return -sigma;
  } else {
    return Dim<1>::SymTensor();
  }
}

inline
Dim<2>::SymTensor
tensileStressCorrection(const Dim<2>::SymTensor& sigma) {
  const EigenStruct<2> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  Dim<2>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,
                           0.0, (lambday > 0.0 ? -lambday : 0.0));
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
  Dim<3>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0, 0.0,
                           0.0, (lambday > 0.0 ? -lambday : 0.0), 0.0,
                           0.0, 0.0, (lambdaz > 0.0 ? -lambdaz : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SolidCRKSPHHydroBaseRZ::
SolidCRKSPHHydroBaseRZ(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                       ArtificialViscosity<Dimension>& Q,
                       const TableKernel<Dimension>& W,
                       const TableKernel<Dimension>& WPi,
                       const double filter,
                       const double cfl,
                       const bool useVelocityMagnitudeForDt,
                       const bool compatibleEnergyEvolution,
                       const bool evolveTotalEnergy,
                       const bool XSPH,
                       const MassDensityType densityUpdate,
                       const HEvolutionType HUpdate,
                       const CRKOrder correctionOrder,
                       const CRKVolumeType volumeType,
                       const double epsTensile,
                       const double nTensile,
                       const bool damageRelieveRubble):
  SolidCRKSPHHydroBase<Dimension>(smoothingScaleMethod, 
                                  Q,
                                  W,
                                  WPi,
                                  filter,
                                  cfl,
                                  useVelocityMagnitudeForDt,
                                  compatibleEnergyEvolution,
                                  evolveTotalEnergy,
                                  XSPH,
                                  densityUpdate,
                                  HUpdate,
                                  correctionOrder,
                                  volumeType,
                                  epsTensile,
                                  nTensile,
                                  damageRelieveRubble),
  mDeviatoricStressTT(FieldStorageType::CopyFields),
  mDdeviatoricStressTTDt(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SolidCRKSPHHydroBaseRZ::
~SolidCRKSPHHydroBaseRZ() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
initializeProblemStartup(DataBase<Dim<2> >& dataBase) {

  dataBase.isRZ = true;

  // Correct the mass to mass/r.
  auto mass = dataBase.fluidMass();
  const auto pos = dataBase.fluidPosition();
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) /= circi;
    }
  }

  // Call the ancestor.
  SolidCRKSPHHydroBase<Dimension>::initializeProblemStartup(dataBase);

  // Convert back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) *= circi;
    }
  }

  // Create storage for the state we're holding.
  mDeviatoricStressTT = dataBase.newSolidFieldList(0.0, SolidFieldNames::deviatoricStressTT);
  mDdeviatoricStressTTDt = dataBase.newSolidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + SolidFieldNames::deviatoricStressTT);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
registerState(DataBase<Dim<2> >& dataBase,
              State<Dim<2> >& state) {

  typedef State<Dimension>::PolicyPointer PolicyPointer;

  // Call the ancestor.
  SolidCRKSPHHydroBase<Dimension>::registerState(dataBase, state);

  // Create the local storage.
  dataBase.resizeFluidFieldList(mDeviatoricStressTT, 0.0, SolidFieldNames::deviatoricStressTT, false);

  // Reregister the deviatoric stress and plastic strain policies to the RZ specialized versions
  // that account for the theta-theta component of the stress.
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> ps = state.fields(SolidFieldNames::plasticStrain, 0.0);
  PolicyPointer deviatoricStressPolicy(new RZDeviatoricStressPolicy());
  PolicyPointer plasticStrainPolicy(new RZPlasticStrainPolicy());
  state.enroll(S, deviatoricStressPolicy);
  state.enroll(ps, plasticStrainPolicy);
  state.enroll(mDeviatoricStressTT);

  // Are we using the compatible energy evolution scheme?
  // If so we need to override the ordinary energy registration with a specialized version.
  if (mCompatibleEnergyEvolution) {
    FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
    PolicyPointer thermalEnergyPolicy(new NonSymmetricSpecificThermalEnergyPolicy<Dimension>(dataBase));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);

    // Get the policy for the position, and add the specific energy as a dependency.
    PolicyPointer positionPolicy = state.policy(state.buildFieldKey(HydroFieldNames::position, UpdatePolicyBase<Dimension>::wildcard()));
    positionPolicy->addDependency(HydroFieldNames::specificThermalEnergy);
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
registerDerivatives(DataBase<Dim<2> >& dataBase,
                    StateDerivatives<Dim<2> >& derivs) {

  // Call the ancestor method.
  SolidCRKSPHHydroBase<Dimension>::registerDerivatives(dataBase, derivs);

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  const string DSTTDtName = IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStressTT;
  dataBase.resizeFluidFieldList(mDdeviatoricStressTTDt, 0.0, DSTTDtName, false);

  derivs.enroll(mDdeviatoricStressTTDt);
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
finalize(const Dim<2>::Scalar time,
         const Dim<2>::Scalar dt,
         DataBase<Dim<2> >& dataBase,
         State<Dim<2> >& state,
         StateDerivatives<Dim<2> >& derivs) {

  // Convert the mass to mass per unit length first.
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) /= circi;
    }
  }

  // Base class finalization does most of the work.
  SolidCRKSPHHydroBase<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Now convert back to true masses.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto& xi = pos(nodeListi, i);
      const auto circi = 2.0*M_PI*abs(xi.y());
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
applyGhostBoundaries(State<Dim<2> >& state,
                     StateDerivatives<Dim<2> >& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary BCs.
  SolidCRKSPHHydroBase<Dim<2> >::applyGhostBoundaries(state, derivs);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Scale back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) *= circi;
    }
  }

  // Apply boundary conditions to our extra strength variables.
  FieldList<Dimension, Scalar> STT = state.fields(SolidFieldNames::deviatoricStressTT, 0.0);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(STT);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
enforceBoundaries(State<Dim<2> >& state,
                  StateDerivatives<Dim<2> >& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary BCs.
  SolidCRKSPHHydroBase<Dim<2> >::enforceBoundaries(state, derivs);

  // Scale back to mass.
  // We also ensure no point approaches the z-axis too closely.
  FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    const Scalar nPerh = mass[nodeListi]->nodeList().nodesPerSmoothingScale();
    for (unsigned i = 0; i != n; ++i) {
      Vector& posi = pos(nodeListi, i);
      const Scalar circi = 2.0*M_PI*abs(posi.y());
      mass(nodeListi, i) *= circi;
    }
  }

  // Enforce boundary conditions on the extra strength variables.
  FieldList<Dimension, Scalar> STT = state.fields(SolidFieldNames::deviatoricStressTT, 0.0);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(STT);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
dumpState(FileIO& file, const string& pathName) const {

  // Ancestor method.
  SolidCRKSPHHydroBase<Dimension>::dumpState(file, pathName);

  file.write(mDeviatoricStressTT, pathName + "/deviatoricStressTT");
  file.write(mDdeviatoricStressTTDt, pathName + "/DdeviatoricStressTTDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
restoreState(const FileIO& file, const string& pathName) {
 
  // Ancestor method.
  SolidCRKSPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mDeviatoricStressTT, pathName + "/deviatoricStressTT");
  file.read(mDdeviatoricStressTTDt, pathName + "/DdeviatoricStressTTDt");
}

}

#ifdef _OPENMP
#include "SolidCRKSPHEvaluateDerivativesRZ_OpenMP.cc"
#else
#include "SolidCRKSPHEvaluateDerivativesRZ.cc"
#endif
