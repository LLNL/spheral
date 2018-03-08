//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBaseRZ -- The axisymmetric (RZ) SPH/ASPH solid material
//                        hydrodynamic package for Spheral++.
//
// This RZ version is a naive area-weighting implementation, nothing as
// highfalutin as the Garcia-Senz approach.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Mon May  9 11:01:51 PDT 2016
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DamagedNodeCouplingWithFrags.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/RZDeviatoricStressPolicy.hh"
#include "Strength/RZPlasticStrainPolicy.hh"
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
#include "SolidMaterial/SolidEquationOfState.hh"

#include "SolidSPHHydroBaseRZ.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

namespace Spheral {
namespace SPHSpace {

using namespace std;
using SPHSpace::SPHHydroBase;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NodeSpace::SolidNodeList;
using SolidMaterial::SolidEquationOfState;
using FileIOSpace::FileIO;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Compute the artificial tensile stress correction tensor for the given 
// stress tensor
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SolidSPHHydroBaseRZ::
SolidSPHHydroBaseRZ(const SmoothingScaleBase<Dim<2> >& smoothingScaleMethod,
                    ArtificialViscosity<Dim<2> >& Q,
                    const TableKernel<Dim<2> >& W,
                    const TableKernel<Dim<2> >& WPi,
                    const TableKernel<Dim<2> >& WGrad,
                    const double filter,
                    const double cfl,
                    const bool useVelocityMagnitudeForDt,
                    const bool compatibleEnergyEvolution,
                    const bool evolveTotalEnergy,
                    const bool gradhCorrection,
                    const bool XSPH,
                    const bool correctVelocityGradient,
                    const bool sumMassDensityOverAllNodeLists,
                    const PhysicsSpace::MassDensityType densityUpdate,
                    const PhysicsSpace::HEvolutionType HUpdate,
                    const double epsTensile,
                    const double nTensile,
                    const Vector& xmin,
                    const Vector& xmax):
  SolidSPHHydroBase<Dim<2> >(smoothingScaleMethod, 
                             Q,
                             W,
                             WPi,
                             WGrad,
                             filter,
                             cfl,
                             useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution,
                             evolveTotalEnergy,
                             gradhCorrection,
                             XSPH,
                             correctVelocityGradient,
                             sumMassDensityOverAllNodeLists,
                             densityUpdate,
                             HUpdate,
                             epsTensile,
                             nTensile,
                             xmin,
                             xmax),
  mDeviatoricStressTT(FieldSpace::FieldStorageType::CopyFields),
  mDdeviatoricStressTTDt(FieldSpace::FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SolidSPHHydroBaseRZ::
~SolidSPHHydroBaseRZ() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
void
SolidSPHHydroBaseRZ::
initializeProblemStartup(DataBase<Dim<2> >& dataBase) {

  // Call the ancestor.
  SolidSPHHydroBase<Dim<2> >::initializeProblemStartup(dataBase);

  dataBase.isRZ = true;

  // Create storage for the state we're holding.
  mDeviatoricStressTT = dataBase.newSolidFieldList(0.0, SolidFieldNames::deviatoricStressTT);
  mDdeviatoricStressTTDt = dataBase.newSolidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + SolidFieldNames::deviatoricStressTT);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SolidSPHHydroBaseRZ::
registerState(DataBase<Dim<2> >& dataBase,
              State<Dim<2> >& state) {

  typedef State<Dimension>::PolicyPointer PolicyPointer;

  // Call the ancestor.
  SolidSPHHydroBase<Dim<2> >::registerState(dataBase, state);

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
SolidSPHHydroBaseRZ::
registerDerivatives(DataBase<Dim<2> >& dataBase,
                    StateDerivatives<Dim<2> >& derivs) {

  // Call the ancestor method.
  SolidSPHHydroBase<Dim<2> >::registerDerivatives(dataBase, derivs);

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
SolidSPHHydroBaseRZ::
finalize(const Dim<2>::Scalar time,
         const Dim<2>::Scalar dt,
         DataBase<Dim<2> >& dataBase,
         State<Dim<2> >& state,
         StateDerivatives<Dim<2> >& derivs) {

  // If we're going to do the SPH summation density, we need to convert the mass
  // to mass per unit length first.
  if (densityUpdate() == PhysicsSpace::MassDensityType::RigorousSumDensity or
      densityUpdate() == PhysicsSpace::MassDensityType::CorrectedSumDensity) {
    FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
    const unsigned numNodeLists = mass.numFields();
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = mass[nodeListi]->numElements();
      for (unsigned i = 0; i != n; ++i) {
        const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
        mass(nodeListi, i) /= circi;
      }
    }
  }

  // Base class finalization does most of the work.
  SPHHydroBase<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Now convert back to true masses and mass densities.  We also apply the RZ
  // correction factor to the mass density.
  if (densityUpdate() == PhysicsSpace::MassDensityType::RigorousSumDensity or
      densityUpdate() == PhysicsSpace::MassDensityType::CorrectedSumDensity) {
    const TableKernel<Dimension>& W = this->kernel();
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const unsigned numNodeLists = massDensity.numFields();
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = massDensity[nodeListi]->numElements();
      for (unsigned i = 0; i != n; ++i) {
        const Vector& xi = position(nodeListi, i);
        const Scalar circi = 2.0*M_PI*abs(xi.y());
        mass(nodeListi, i) *= circi;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidSPHHydroBaseRZ::
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

  // Apply ordinary SPH BCs.
  SolidSPHHydroBase<Dim<2> >::applyGhostBoundaries(state, derivs);
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
SolidSPHHydroBaseRZ::
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

  // Apply ordinary SPH BCs.
  SolidSPHHydroBase<Dim<2> >::enforceBoundaries(state, derivs);

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
SolidSPHHydroBaseRZ::
dumpState(FileIO& file, const string& pathName) const {

  // Ancestor method.
  SolidSPHHydroBase<Dimension>::dumpState(file, pathName);

  file.write(mDeviatoricStressTT, pathName + "/deviatoricStressTT");
  file.write(mDdeviatoricStressTTDt, pathName + "/DdeviatoricStressTTDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
void
SolidSPHHydroBaseRZ::
restoreState(const FileIO& file, const string& pathName) {
 
  // Ancestor method.
  SolidSPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mDeviatoricStressTT, pathName + "/deviatoricStressTT");
  file.read(mDdeviatoricStressTTDt, pathName + "/DdeviatoricStressTTDt");
}

}
}

// Include the appropriate evaluateDerivatvies.
#ifdef _OPENMP
#include "SolidSPHEvaluateDerivativesRZ_OpenMP.cc"
#else
#include "SolidSPHEvaluateDerivativesRZ.cc"
#endif
