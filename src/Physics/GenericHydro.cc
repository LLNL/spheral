//---------------------------------Spheral++----------------------------------//
// GenericHydro -- The base class for all Spheral++ hydro implementations.
//
// Created by JMO, Sat May 20 22:50:20 PDT 2000
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>

#include <sstream>
#include <algorithm>

using namespace std;

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

//------------------------------------------------------------------------------
// Helper method to compute the magnitude of the shear term in the given
// dvdx tensor.
//------------------------------------------------------------------------------
double
computeShearMagnitude(const Dim<1>::Tensor& dvdx) {
  return 0.0;
}

double
computeShearMagnitude(const Dim<2>::Tensor& dvdx) {
  return std::abs(dvdx(1,0) - dvdx(0,1));
}

double
computeShearMagnitude(const Dim<3>::Tensor& dvdx) {
  return sqrt(FastMath::square(dvdx(2,1) - dvdx(1,2)) +
              FastMath::square(dvdx(2,0) - dvdx(0,2)) +
              FastMath::square(dvdx(1,0) - dvdx(0,1)));
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

  // Get some useful fluid variables from the DataBase.
  const FieldList<Dimension, int> mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  const FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity());
  const int numNodeLists = connectivityMap.nodeLists().size();

  // Initialize the return value to some impossibly high value.
  Scalar minDt = FLT_MAX;

  // Set up some history variables to track what set's our minimum Dt.
  Scalar lastMinDt = minDt;
  int lastNodeID;
  string lastNodeListName, reason;
  Scalar lastNodeScale, lastCs, lastRho, lastEps, lastDivVelocity, lastShearVelocity;
  Vector lastVelocity, lastAcc;

  // Loop over every fluid node.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
       nodeListItr != dataBase.fluidNodeListEnd();
       ++nodeListItr, ++nodeListi) {
    const FluidNodeList<Dimension>& fluidNodeList = **nodeListItr;
    const Scalar nPerh = fluidNodeList.nodesPerSmoothingScale();
    CHECK(nPerh > 0.0);
    // const Scalar kernelExtent = fluidNodeList.neighbor().kernelExtent();
    // CHECK(kernelExtent > 0.0);

    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // If this node is masked, don't worry about it.
      if (mask(nodeListi, i) == 1) {

        // Get this nodes minimum characteristic smoothing scale.
        CHECK2(H(nodeListi, i).Determinant() >  0.0,
               "Bad H tensor : " << H(nodeListi, i) << " : " << fluidNodeList.name() << " " << i << " " << fluidNodeList.firstGhostNode());
        const Scalar nodeScale = 1.0/H(nodeListi, i).eigenValues().maxElement()/nPerh;
        // const Scalar nodeScale = 1.0/Dimension::rootnu(H(nodeListi, i).Determinant());
        //     const Scalar nodeScale = nodeExtent(nodeListi, i).minElement()/kernelExtent;

        // Sound speed limit.
        const double csDt = nodeScale/(soundSpeed(nodeListi, i) + FLT_MIN);
        if (csDt < minDt) {
          minDt = csDt;
          reason = "sound speed limit";
        }

        // Artificial viscosity effective sound speed.
        CHECK(rho(nodeListi, i) > 0.0);
        const Scalar csq = sqrt(maxViscousPressure(nodeListi, i)/rho(nodeListi, i));
        const double csqDt = nodeScale/(csq + FLT_MIN);
        if (csqDt < minDt) {
          minDt = csqDt;
          reason = "artificial viscosity sound speed limit";
        }

        // Velocity divergence limit.
        const Scalar divVelocity = DvDx(nodeListi, i).Trace();
        const double divvDt = 1.0/(std::abs(divVelocity) + FLT_MIN);
        if (divvDt < minDt) {
          minDt = divvDt;
          reason = "velocity divergence";
        }

        // Maximum velocity difference limit.
        const Vector& vi = velocity(nodeListi, i);
        const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          const unsigned firstGhostNodej = velocity[nodeListj]->nodeList().firstGhostNode();
          const vector<int>& connectivity = fullConnectivity[nodeListj];
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            const Vector& vj = velocity(nodeListj, j);
            const Scalar vij = (vi - vj).magnitude();
            const Scalar dtVelDiff = nodeScale/max(1.0e-100, vij);
            if (dtVelDiff < minDt) {
              minDt = dtVelDiff;
              reason = "pairwise velocity difference";
            }
          }
        }

        //     // Eigenvalues of the stress-strain tensor.
        //     Vector eigenValues = DvDx(nodeListi, i).Symmetric().eigenValues();
        //     Scalar maxValue = -1.0;
        //     for (int i = 0; i < Dimension::nDim; ++i) {
        //       maxValue = max(maxValue, fabs(eigenValues(i)));
        //     }
        //     const double strainDt = 1.0/(maxValue + FLT_MIN);
        //     if (strainDt < minDt) {
        //       minDt = strainDt;
        //       reason = "strain limit";
        //     }

        //     // Limit by the velocity shear.
        //     const Scalar shearVelocity = computeShearMagnitude(DvDx(nodeListi, i));
        //     const double shearDt = 1.0/(shearVelocity + FLT_MIN);
        //     if (shearDt < minDt) {
        //       minDt = shearDt;
        //       reason = "velocity shear limit";
        //     }

        // Total acceleration limit.
        const double dtAcc = sqrt(nodeScale/(DvDt(nodeListi, i).magnitude() + FLT_MIN));
        if (dtAcc < minDt) {
          minDt = dtAcc;
          reason = "total acceleration";
        }

        // If requested, limit against the absolute velocity.
        if (useVelocityMagnitudeForDt()) {
          const double velDt = nodeScale/(velocity(nodeListi, i).magnitude() + 1.0e-10);
          if (velDt < minDt) {
            minDt = velDt;
            reason = "velocity magnitude";
          }
        }

        if (minDt < lastMinDt) {
          lastMinDt = minDt;
          lastNodeID = i;
          lastNodeListName = fluidNodeList.name();
          lastNodeScale = nodeScale;
          lastCs = soundSpeed(nodeListi, i);
          lastAcc = DvDt(nodeListi, i);
          //      lastCsq = csq;
          lastRho = rho(nodeListi, i);
          lastEps = eps(nodeListi, i);
          lastVelocity = velocity(nodeListi, i);
          lastDivVelocity = divVelocity;
          //       lastShearVelocity = shearVelocity;
        }
      }
    }
  }

  stringstream reasonStream;
  reasonStream << "mindt = " << minDt << endl
	       << reason << endl
               << "  (nodeList, node) = (" << lastNodeListName << ", " << lastNodeID << ") | "
               << "  h = " << lastNodeScale << endl
               << "  cs = " << lastCs << endl
               << "  acc = " << lastAcc << endl
//               << " csq = " << lastCsq << endl
               << "  rho = " << lastRho << endl
               << "  eps = " << lastEps << endl
               << "  velocity = " << lastVelocity << endl
               << "  dtcs = " << lastNodeScale/(lastCs + FLT_MIN) << endl
               << "  divVelocity = " << lastDivVelocity << endl
//               << "  dtcsq = " << lastNodeScale/(lastCsq + FLT_MIN) << endl
//                << "  dtdivV = " << 1.0/(fabs(lastDivVelocity) + FLT_MIN) << endl
//                << "  shearVelocity = " << lastShearVelocity << endl
               << ends;

  // Now build the result.
  TimeStepType result(cfl()*minDt, reasonStream.str());

  return result;
}

}
}

