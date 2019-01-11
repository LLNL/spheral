//---------------------------------Spheral++----------------------------------//
// SPHHydroBaseRZ -- An SPH/ASPH hydrodynamic package for Spheral++,
//                   specialized for 2D RZ (cylindrical) geometry.
//
// This RZ version is a naive area-weighting implementation, nothing as
// highfalutin as the Garcia-Senz approach.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Fri May  6 16:18:36 PDT 2016
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computeSPHOmegaGradhCorrection.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"

#include "SPHHydroBaseRZ.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
using std::vector;
using std::map;
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
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SPHHydroBaseRZ::
SPHHydroBaseRZ(const SmoothingScaleBase<Dim<2> >& smoothingScaleMethod,
               ArtificialViscosity<Dim<2> >& Q,
               const TableKernel<Dim<2> >& W,
               const TableKernel<Dim<2> >& WPi,
               const double filter,
               const double cfl,
               const bool useVelocityMagnitudeForDt,
               const bool compatibleEnergyEvolution,
               const bool evolveTotalEnergy,
               const bool gradhCorrection,
               const bool XSPH,
               const bool correctVelocityGradient,
               const bool sumMassDensityOverAllNodeLists,
               const MassDensityType densityUpdate,
               const HEvolutionType HUpdate,
               const double epsTensile,
               const double nTensile,
               const Vector& xmin,
               const Vector& xmax):
  SPHHydroBase<Dim<2> >(smoothingScaleMethod,
                        Q,
                        W,
                        WPi,
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
                        xmax) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SPHHydroBaseRZ::
~SPHHydroBaseRZ() {
}

//------------------------------------------------------------------------------
// On problem start up, we set the RZ flag on the DataBase.
//------------------------------------------------------------------------------
void
SPHHydroBaseRZ::
initializeProblemStartup(DataBase<Dim<2> >& dataBase) {
  dataBase.isRZ = true;
  SPHHydroBase<Dim<2> >::initializeProblemStartup(dataBase);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SPHHydroBaseRZ::
registerState(DataBase<Dim<2> >& dataBase,
              State<Dim<2> >& state) {

  typedef State<Dimension>::PolicyPointer PolicyPointer;

  // The base class does most of it.
  SPHHydroBase<Dim<2> >::registerState(dataBase, state);

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
// Finalize the hydro.
//------------------------------------------------------------------------------
void
SPHHydroBaseRZ::
finalize(const Dim<2>::Scalar time,
         const Dim<2>::Scalar dt,
         DataBase<Dim<2> >& dataBase,
         State<Dim<2> >& state,
         StateDerivatives<Dim<2> >& derivs) {

  // If we're going to do the SPH summation density, we need to convert the mass
  // to mass per unit length first.
  if (densityUpdate() == MassDensityType::RigorousSumDensity or
      densityUpdate() == MassDensityType::CorrectedSumDensity) {
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
  if (densityUpdate() == MassDensityType::RigorousSumDensity or
      densityUpdate() == MassDensityType::CorrectedSumDensity) {
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
SPHHydroBaseRZ::
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
  SPHHydroBase<Dim<2> >::applyGhostBoundaries(state, derivs);
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
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SPHHydroBaseRZ::
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
  SPHHydroBase<Dim<2> >::enforceBoundaries(state, derivs);

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
}

}

// Include the appropriate evaluateDerivatvies.
#ifdef _OPENMP
#include "SPHEvaluateDerivativesRZ_OpenMP.cc"
#else
#include "SPHEvaluateDerivativesRZ.cc"
#endif
