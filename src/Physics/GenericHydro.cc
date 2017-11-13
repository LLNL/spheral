//---------------------------------Spheral++----------------------------------//
// GenericHydro -- The base class for all Spheral++ hydro implementations.
//
// Created by JMO, Sat May 20 22:50:20 PDT 2000
//----------------------------------------------------------------------------//
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "GenericHydro.hh"

#include "Physics.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Strength/SolidFieldNames.hh"

#include <limits>
#include <algorithm>

using namespace std;

namespace Spheral {
namespace PhysicsSpace {

using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using NeighborSpace::ConnectivityMap;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;

namespace {

//------------------------------------------------------------------------------
// Helper method to compute the magnitude of the shear term in the given
// dvdx tensor.
//------------------------------------------------------------------------------
inline
double
computeShearMagnitude(const Dim<1>::Tensor& dvdx) {
  return 0.0;
}

inline
double
computeShearMagnitude(const Dim<2>::Tensor& dvdx) {
  return std::abs(dvdx(1,0) - dvdx(0,1));
}

inline
double
computeShearMagnitude(const Dim<3>::Tensor& dvdx) {
  return sqrt(FastMath::square(dvdx(2,1) - dvdx(1,2)) +
              FastMath::square(dvdx(2,0) - dvdx(0,2)) +
              FastMath::square(dvdx(1,0) - dvdx(0,1)));
}

}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
GenericHydro<Dimension>::
GenericHydro(const TableKernel<Dimension>& W,
             const TableKernel<Dimension>& WPi,
             ArtificialViscosity<Dimension>& Q,
             const double cfl,
             const bool useVelocityMagnitudeForDt):
  Physics<Dimension>(),
  mArtificialViscosity(Q),
  mKernel(W),
  mPiKernel(WPi),
  mCfl(cfl),
  mUseVelocityMagnitudeForDt(useVelocityMagnitudeForDt),
  mMinMasterNeighbor(INT_MAX),
  mMaxMasterNeighbor(0),
  mSumMasterNeighbor(0),
  mNormMasterNeighbor(0),
  mMinCoarseNeighbor(INT_MAX),
  mMaxCoarseNeighbor(0),
  mSumCoarseNeighbor(0),
  mNormCoarseNeighbor(0),
  mMinRefineNeighbor(INT_MAX),
  mMaxRefineNeighbor(0),
  mSumRefineNeighbor(0),
  mNormRefineNeighbor(0),
  mMinActualNeighbor(INT_MAX),
  mMaxActualNeighbor(0),
  mSumActualNeighbor(0),
  mNormActualNeighbor(0) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GenericHydro<Dimension>::~GenericHydro() {
}

//------------------------------------------------------------------------------
// Determine the timestep requirements for a hydro step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename GenericHydro<Dimension>::TimeStepType
GenericHydro<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   typename Dimension::Scalar currentTime) const {

  const double tiny = numeric_limits<double>::epsilon();

  // Get some useful fluid variables from the DataBase.
  const auto  mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  eps = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  const auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  const auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity());
  const auto  numNodeLists = connectivityMap.nodeLists().size();

  // Check for deviatoric stress.
  const auto haveDS = state.fieldNameRegistered(SolidFieldNames::deviatoricStress);
  FieldList<Dimension, SymTensor> S;
  if (haveDS) S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);

  // Check if the longitudinal sound speed is registered.
  const auto haveLongCs = state.fieldNameRegistered(SolidFieldNames::longitudinalSoundSpeed);
  FieldList<Dimension, Scalar> csl;
  if (haveLongCs) csl = state.fields(SolidFieldNames::longitudinalSoundSpeed, 0.0);

  // Initialize the return value to some impossibly high value.
  auto minDt = make_pair(numeric_limits<double>::max(), string());

  // Loop over every fluid node.
  // #pragma omp declare reduction (MINPAIR : pair<double,string> : omp_out = (omp_out.first < omp_in.first ? omp_out : omp_in)) initializer(omp_priv = pair<double,string>(numeric_limits<double>::max(), string("null")))
  // #pragma omp parallel for reduction(MINPAIR:minDt) collapse(2)
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto& fluidNodeList = **(dataBase.fluidNodeListBegin() + nodeListi);
    const auto nPerh = fluidNodeList.nodesPerSmoothingScale();
    CHECK(nPerh > 0.0);

    // Check if we have a longitudinal sound speed for this material.
    const bool useCsl = haveLongCs and csl.haveNodeList(fluidNodeList);
    const Field<Dimension, Scalar>* cslptr;
    if (useCsl) cslptr = *csl.fieldForNodeList(fluidNodeList);

    // Check if we have a deviatoric stress for this material.
    const bool useS = haveDS and S.haveNodeList(fluidNodeList);
    const Field<Dimension, SymTensor>* Sptr;
    if (useS) Sptr = *S.fieldForNodeList(fluidNodeList);

    // Walk all the nodes in this FluidNodeList.
    const auto ni = connectivityMap.numNodes(nodeListi);
    // #pragma omp parallel for reduction(MINPAIR:minDt)
#pragma omp parallel
    {
      auto minDt_local = minDt;
#pragma omp for
      for (auto k = 0; k < ni; ++k) {
        const auto i = connectivityMap.ithNode(nodeListi, k);

        // If this node is masked, don't worry about it.
        if (mask(nodeListi, i) == 1) {

          // Get this nodes minimum characteristic smoothing scale.
          CHECK2(H(nodeListi, i).Determinant() >  0.0,
                 "Bad H tensor : " << H(nodeListi, i) << " : " << fluidNodeList.name() << " " << i << " " << fluidNodeList.firstGhostNode());
          const Scalar nodeScale = 1.0/H(nodeListi, i).eigenValues().maxElement()/nPerh;
          // const Scalar nodeScale = 1.0/Dimension::rootnu(H(nodeListi, i).Determinant());
          //     const Scalar nodeScale = nodeExtent(nodeListi, i).minElement()/kernelExtent;

          // Sound speed limit.
          const auto csDt = nodeScale/(cs(nodeListi, i) + tiny);
          if (csDt < minDt_local.first) {
            minDt_local = make_pair(csDt, ("Sound speed limit: dt = " + to_string(csDt) + "\n" +
                                     "                   cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                     "            nodeScale = " + to_string(nodeScale)));
          }

          // Longitudinal sound speed limit.
          if (useCsl) {
            const auto csDt = nodeScale/((*cslptr)(i) + tiny);
            if (csDt < minDt_local.first) {
              minDt_local = make_pair(csDt, ("Longitudinal sound speed limit: dt = " + to_string(csDt) + "\n" + 
                                       "                                cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                       "                               csl = " + to_string((*cslptr)(i)) + "\n" +
                                       "                         nodeScale = " + to_string(nodeScale)));
            }
          }

          // Deviatoric stress limit.
          if (useS) {
            const auto csS = sqrt((*Sptr)(i).eigenValues().maxAbsElement()/rho(nodeListi, i));
            const auto SDt = nodeScale/(csS + tiny);
            if (SDt < minDt_local.first) {
              minDt_local = make_pair(SDt, ("Deviatoric stress effective sound speed limit: dt = " + to_string(SDt) + "\n" +
                                      "                                               cs = " + to_string(cs(nodeListi, i)) + "\n" + 
                                      "                                              csS = " + to_string(csS) + "\n" +
                                      "                                              rho = " + to_string(rho(nodeListi, i)) + "\n" +
                                      "                                        nodeScale = " + to_string(nodeScale)));
            }
          }

          // Artificial viscosity effective sound speed.
          CHECK(rho(nodeListi, i) > 0.0);
          const auto csq = sqrt(maxViscousPressure(nodeListi, i)/rho(nodeListi, i));
          const auto csqDt = nodeScale/(csq + tiny);
          if (csqDt < minDt_local.first) {
            minDt_local = make_pair(csqDt, ("Artificial viscosity sound speed limit: dt = " + to_string(csqDt) + "\n" + 
                                      "                                        cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                      "                                       csQ = " + to_string(csq) + "\n" +
                                      "                                       rho = " + to_string(rho(nodeListi, i)) + "\n" +
                                      "                                 nodeScale = " + to_string(nodeScale)));
          }

          // Velocity divergence limit.
          const auto divVelocity = DvDx(nodeListi, i).Trace();
          const auto divvDt = 1.0/(std::abs(divVelocity) + tiny);
          if (divvDt < minDt_local.first) {
            minDt_local = make_pair(divvDt, ("Velocity divergence limit: dt = " + to_string(divvDt) + "\n" +
                                       "                 div velocity = " + to_string(divVelocity)));
          }

          // Maximum velocity difference limit.
          const auto& xi = position(nodeListi, i);
          const auto& vi = velocity(nodeListi, i);
          const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
          for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
            const auto& connectivity = fullConnectivity[nodeListj];
            for (auto jItr = connectivity.begin();
                 jItr != connectivity.end();
                 ++jItr) {
              const auto  j = *jItr;
              const auto& xj = position(nodeListj, j);
              const auto& vj = velocity(nodeListj, j);
              // const Scalar vij = std::abs((vj - vi).dot((xi - xj).unitVector()));
              const auto  vij = (vi - vj).magnitude();
              const auto  dtVelDiff = nodeScale*safeInvVar(vij, 1e-30);
              if (dtVelDiff < minDt_local.first) {
                minDt_local = make_pair(dtVelDiff, ("Pairwise velocity difference limit: dt = " + to_string(dtVelDiff) + "\n" + 
                                              "                        (nodeListi, i) = " + to_string(nodeListi) + " " + to_string(i) + "\n" +
                                              "                        (nodeListj, j) = " + to_string(nodeListj) + " " + to_string(j) + "\n" +
                                              "                                   vij = " + to_string(vij) + "\n" +
                                              "                             nodeScale = " + to_string(nodeScale)));
              }
            }
          }

          //     // Eigenvalues of the stress-strain tensor.
          //     Vector eigenValues = DvDx(nodeListi, i).Symmetric().eigenValues();
          //     Scalar maxValue = -1.0;
          //     for (int i = 0; i < Dimension::nDim; ++i) {
          //       maxValue = max(maxValue, fabs(eigenValues(i)));
          //     }
          //     const double strainDt = 1.0/(maxValue + tiny);
          //     if (strainDt < minDt) {
          //       minDt = strainDt;
          //       reason = "strain limit";
          //     }

          //     // Limit by the velocity shear.
          //     const Scalar shearVelocity = computeShearMagnitude(DvDx(nodeListi, i));
          //     const double shearDt = 1.0/(shearVelocity + tiny);
          //     if (shearDt < minDt) {
          //       minDt = shearDt;
          //       reason = "velocity shear limit";
          //     }

          // Total acceleration limit.
          const auto dtAcc = sqrt(nodeScale/(DvDt(nodeListi, i).magnitude() + tiny));
          if (dtAcc < minDt_local.first) {
            minDt_local = make_pair(dtAcc, ("Total acceleration limit: dt = " + to_string(dtAcc) + "\n" + 
                                      "              |acceleration| = " + to_string(DvDt(nodeListi, i).magnitude()) + "\n" +
                                      "                   nodeScale = " + to_string(nodeScale)));
          }

          // If requested, limit against the absolute velocity.
          if (useVelocityMagnitudeForDt()) {
            const auto velDt = nodeScale/(velocity(nodeListi, i).magnitude() + 1.0e-10);
            if (velDt < minDt_local.first) {
              minDt_local = make_pair(velDt, ("Velocity magnitude limit: dt = " + to_string(velDt) + "\n" +
                                        "                        |vi| = " + to_string(velocity(nodeListi, i).magnitude()) + "\n" +
                                        "                   nodeScale = " + to_string(nodeScale)));
            }
          }
        }
      }

#pragma omp critical
      if (minDt_local.first < minDt.first) minDt = minDt_local;
    }
  }

  // Scale by the cfl safety factor.
  minDt.first *= cfl();
  return minDt;
}

}
}

