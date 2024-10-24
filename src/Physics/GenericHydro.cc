//---------------------------------Spheral++----------------------------------//
// GenericHydro -- The base class for all Spheral++ hydro implementations.
//
// Created by JMO, Sat May 20 22:50:20 PDT 2000
//----------------------------------------------------------------------------//
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "GenericHydro.hh"

#include "Physics.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Strength/SolidFieldNames.hh"
#include "Utilities/range.hh"

#include <limits>
#include <algorithm>
#include <sstream>

using std::vector;
using std::string;
using std::to_string;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::make_pair;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Helper method to compute the magnitude of the shear term in the given
// dvdx tensor.
//------------------------------------------------------------------------------
inline
double
computeShearMagnitude(const Dim<1>::Tensor& /*dvdx*/) {
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

// Provide to_string for Spheral Vector
template<typename Vector>
std::string
vec_to_string(const Vector& vec) {
  std::ostringstream oss;
  oss << vec << std::endl;
  return oss.str();
}

}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
GenericHydro<Dimension>::
GenericHydro(ArtificialViscosity<Dimension>& Q,
             const double cfl,
             const bool useVelocityMagnitudeForDt):
  Physics<Dimension>(),
  mArtificialViscosity(Q),
  mCfl(cfl),
  mUseVelocityMagnitudeForDt(useVelocityMagnitudeForDt),
  mMinMasterNeighbor(INT_MAX),
  mMaxMasterNeighbor(0),
  mSumMasterNeighbor(0),
  mMinCoarseNeighbor(INT_MAX),
  mMaxCoarseNeighbor(0),
  mSumCoarseNeighbor(0),
  mMinRefineNeighbor(INT_MAX),
  mMaxRefineNeighbor(0),
  mSumRefineNeighbor(0),
  mMinActualNeighbor(INT_MAX),
  mMaxActualNeighbor(0),
  mSumActualNeighbor(0),
  mNormMasterNeighbor(0),
  mNormCoarseNeighbor(0),
  mNormRefineNeighbor(0),
  mNormActualNeighbor(0),
  mDTrank(0u),
  mDTNodeList(0u),
  mDTnode(0u),
  mDTreason() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GenericHydro<Dimension>::~GenericHydro() {
}

//------------------------------------------------------------------------------
// Determine the timestep requirements for a hydro step (explicit)
//------------------------------------------------------------------------------
template<typename Dimension>
typename GenericHydro<Dimension>::TimeStepType
GenericHydro<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   typename Dimension::Scalar /*currentTime*/) const {

  const double tiny = std::numeric_limits<double>::epsilon();

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
  const auto& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity(),
                                                         this->requireOverlapConnectivity(),
                                                         this->requireIntersectionConnectivity());
  const auto  pairs = connectivityMap.nodePairList();

  // Check for deviatoric stress.
  const auto haveDS = state.fieldNameRegistered(SolidFieldNames::deviatoricStress);
  FieldList<Dimension, SymTensor> S;
  if (haveDS) S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);

  // Check if the longitudinal sound speed is registered.
  const auto haveLongCs = state.fieldNameRegistered(SolidFieldNames::longitudinalSoundSpeed);
  FieldList<Dimension, Scalar> csl;
  if (haveLongCs) csl = state.fields(SolidFieldNames::longitudinalSoundSpeed, 0.0);

  // Initialize the return value to some impossibly high value.
  auto minDt = make_pair(std::numeric_limits<double>::max(), string());

  // Define a function for computing the velocity divergence (different for curvilinear coordinates)
  const auto Fdiv = (GeometryRegistrar::coords() == CoordinateType::Spherical ? +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi[0] + 2.0*veli[0]*safeInv(posi[0]); } :
                     GeometryRegistrar::coords() == CoordinateType::RZ        ? +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi.Trace() + veli[1]*safeInv(posi[1]); } :
                                                                                +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi.Trace(); });
  // Loop over every fluid node.
  // #pragma omp declare reduction (MINPAIR : pair<double,string> : omp_out = (omp_out.first < omp_in.first ? omp_out : omp_in)) initializer(omp_priv = pair<double,string>(std::numeric_limits<double>::max(), string("null")))
  // #pragma omp parallel for reduction(MINPAIR:minDt) collapse(2)
  for (auto [nodeListi, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListPtrs())) {
    const auto nPerh = fluidNodeListPtr->nodesPerSmoothingScale();
    CHECK(nPerh > 0.0);

    // Check if we have a longitudinal sound speed for this material.
    const auto useCsl = haveLongCs and csl.haveNodeList(*fluidNodeListPtr);
    const Field<Dimension, Scalar>* cslptr = nullptr;
    if (useCsl) cslptr = *csl.fieldForNodeList(*fluidNodeListPtr);

    // Check if we have a deviatoric stress for this material.
    const auto useS = haveDS and S.haveNodeList(*fluidNodeListPtr);
    const Field<Dimension, SymTensor>* Sptr = nullptr;
    if (useS) Sptr = *S.fieldForNodeList(*fluidNodeListPtr);

    // Walk all the nodes in this FluidNodeList.
    const auto ni = connectivityMap.numNodes(nodeListi);
    const auto rank = Process::getRank();
    // #pragma omp parallel for reduction(MINPAIR:minDt)
#pragma omp parallel
    {
      auto minDt_local = minDt;
      auto DTrank_local = mDTrank;
      auto DTNodeList_local = mDTNodeList;
      auto DTnode_local = mDTnode;
      auto DTreason_local = mDTreason;
#pragma omp for
      for (auto k = 0u; k < ni; ++k) {
        const auto i = connectivityMap.ithNode(nodeListi, k);

        // If this node is masked, don't worry about it.
        if (mask(nodeListi, i) == 1) {

          // Get this nodes minimum characteristic smoothing scale.
          const auto& vi = velocity(nodeListi, i);
          const auto& Hi = H(nodeListi, i);
          CHECK2(Hi.Determinant() >  0.0,
                 "Bad H tensor : " << Hi << " : " << fluidNodeListPtr->name() << " " << i << " " << fluidNodeListPtr->firstGhostNode());
          const Scalar nodeScalei = 1.0/Hi.eigenValues().maxElement()/nPerh;
          // const Scalar nodeScalei = 1.0/Dimension::rootnu(Hi.Determinant());
          //     const Scalar nodeScalei = nodeExtent(nodeListi, i).minElement()/kernelExtent;

          // Sound speed limit.
          const auto csi = cs(nodeListi, i);
          const auto csDt = nodeScalei/(csi + tiny);
          if (csDt < minDt_local.first) {
            minDt_local = make_pair(csDt, ("Sound speed limit: dt = " + to_string(csDt) + "\n" +
                                           "                   cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                           "            nodeScale = " + to_string(nodeScalei) + "\n" +
                                           "             material = " + fluidNodeListPtr->name() + "\n" +
                                           "(nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                           "           @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "sound speed";
          }

          // Longitudinal sound speed limit.
          if (useCsl) {
            const auto csDt = nodeScalei/((*cslptr)(i) + tiny);
            if (csDt < minDt_local.first) {
              minDt_local = make_pair(csDt, ("Longitudinal sound speed limit: dt = " + to_string(csDt) + "\n" + 
                                             "                                cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                             "                               csl = " + to_string((*cslptr)(i)) + "\n" +
                                             "                         nodeScale = " + to_string(nodeScalei) + "\n" +
                                             "                          material = " + fluidNodeListPtr->name() + "\n" +
                                             "             (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                             "                        @ position = " + vec_to_string(position(nodeListi, i))));
              DTrank_local = rank;
              DTNodeList_local = nodeListi;
              DTnode_local = i;
              DTreason_local = "longitudinal sound speed";
            }
          }

          // Deviatoric stress limit.
          const auto rhoi = rho(nodeListi, i);
          CHECK(rhoi > 0.0);
          if (useS) {
            const auto csS = sqrt((*Sptr)(i).eigenValues().maxAbsElement()/rhoi);
            const auto SDt = nodeScalei/(csS + tiny);
            if (SDt < minDt_local.first) {
              minDt_local = make_pair(SDt, ("Deviatoric stress effective sound speed limit: dt = " + to_string(SDt) + "\n" +
                                            "                                               cs = " + to_string(cs(nodeListi, i)) + "\n" + 
                                            "                                              csS = " + to_string(csS) + "\n" +
                                            "                                              rho = " + to_string(rhoi) + "\n" +
                                            "                                        nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                                         material = " + fluidNodeListPtr->name() + "\n" +
                                            "                            (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                            "                                       @ position = " + vec_to_string(position(nodeListi, i))));
              DTrank_local = rank;
              DTNodeList_local = nodeListi;
              DTnode_local = i;
              DTreason_local = "deviatoric stress effective sound speed";
            }
          }

          // Artificial viscosity effective sound speed.
          const auto csq = sqrt(maxViscousPressure(nodeListi, i)/rhoi);
          const auto csqDt = nodeScalei/(csq + tiny);
          if (csqDt < minDt_local.first) {
            minDt_local = make_pair(csqDt, ("Artificial viscosity sound speed limit: dt = " + to_string(csqDt) + "\n" + 
                                            "                                        cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                            "                                       csQ = " + to_string(csq) + "\n" +
                                            "                                       rho = " + to_string(rhoi) + "\n" +
                                            "                                 nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                                  material = " + fluidNodeListPtr->name() + "\n" +
                                            "                     (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                            "                                @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "artificial viscosity";
          }

          // Velocity divergence limit.
          const auto divVelocity = Fdiv(DvDx(nodeListi, i), position(nodeListi, i), velocity(nodeListi, i));
          const auto divvDt = 1.0/(std::abs(divVelocity) + tiny);
          if (divvDt < minDt_local.first) {
            minDt_local = make_pair(divvDt, ("Velocity divergence limit: dt = " + to_string(divvDt) + "\n" +
                                             "                 div velocity = " + to_string(divVelocity) + "\n" +
                                             "                     material = " + fluidNodeListPtr->name() + "\n" +
                                             "        (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                             "                   @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "velocity divergence";
          }

          // Total acceleration limit.
          const auto vmagi = vi.magnitude();
          const auto dtAcc = 0.1*std::max(nodeScalei/(vmagi + tiny), vmagi/(DvDt(nodeListi, i).magnitude() + tiny));
          if (dtAcc < minDt_local.first) {
            minDt_local = make_pair(dtAcc, ("Total acceleration limit: dt = " + to_string(dtAcc) + "\n" + 
                                            "              |acceleration| = " + to_string(DvDt(nodeListi, i).magnitude()) + "\n" +
                                            "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                    material = " + fluidNodeListPtr->name() + "\n" +
                                            "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                            "                  @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "acceleration";
          }

          // If requested, limit against the absolute velocity.
          if (useVelocityMagnitudeForDt()) {
            const auto velDt = nodeScalei/(velocity(nodeListi, i).magnitude() + 1.0e-10);
            if (velDt < minDt_local.first) {
              minDt_local = make_pair(velDt, ("Velocity magnitude limit: dt = " + to_string(velDt) + "\n" +
                                              "                        |vi| = " + to_string(velocity(nodeListi, i).magnitude()) + "\n" +
                                              "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                              "                    material = " + fluidNodeListPtr->name() + "\n" +
                                              "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                              "                  @ position = " + vec_to_string(position(nodeListi, i))));
              DTrank_local = rank;
              DTNodeList_local = nodeListi;
              DTnode_local = i;
              DTreason_local = "velocity magnitude";
            }
          }
        }
      }

      // Maximum pair-wise velocity difference limit
      const auto npairs = pairs.size();
#pragma omp for
      for (auto k = 0u; k < npairs; ++k) {
        const auto i = pairs[k].i_node;
        const auto j = pairs[k].j_node;
        const auto nodeListi = pairs[k].i_list;
        const auto nodeListj = pairs[k].j_list;

        const auto& vi = velocity(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  nodeScalei = 1.0/Hi.eigenValues().maxElement()/nPerh;

        const auto& vj = velocity(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  nodeScalej = 1.0/Hj.eigenValues().maxElement()/nPerh;

        const auto  vij = vi - vj;
        const auto  dtVelDiff = std::min(nodeScalei, nodeScalej)*safeInvVar(vij.magnitude(), tiny);
        if (dtVelDiff < minDt_local.first) {
          minDt_local = make_pair(dtVelDiff, ("Pairwise velocity difference limit: dt = " + to_string(dtVelDiff) + "\n" + 
                                              "                              material = " + fluidNodeListPtr->name() + "\n" +
                                              "                  (nodeListi, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                              "                  (nodeListj, i, rank) = (" + to_string(nodeListj) + " " + to_string(j) + " " + to_string(rank) + ")\n" +
                                              "                   @ pos(nodeListi, i) = " + vec_to_string(position(nodeListi, i)) + "\n" +
                                              "                   @ pos(nodeListj, j) = " + vec_to_string(position(nodeListj, j)) + "\n" +
                                              "                                   vij = " + to_string(vij.magnitude()) + "\n" +
                                              "                            nodeScalei = " + to_string(nodeScalei) + "\n" +
                                              "                            nodeScalej = " + to_string(nodeScalej)));
          DTrank_local = rank;
          DTNodeList_local = nodeListi;
          DTnode_local = i;
          DTreason_local = "pairwise velocity difference";
        }
      }

#pragma omp critical
      {
        if (minDt_local.first < minDt.first) {
          minDt = minDt_local;
          mDTrank = DTrank_local;
          mDTNodeList = DTNodeList_local;
          mDTnode = DTnode_local;
          mDTreason = DTreason_local;
        }
      }
    }
  }

  // Scale by the cfl safety factor.
  minDt.first *= cfl();
  return minDt;
}

//------------------------------------------------------------------------------
// Determine the timestep requirements for a hydro step (implicit)
//
// In the implicit case we're going to limit solely by motion, such that
// connectivity/topology shouldn't change too much timestep to timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
typename GenericHydro<Dimension>::TimeStepType
GenericHydro<Dimension>::
dtImplicit(const DataBase<Dimension>& dataBase,
           const State<Dimension>& state,
           const StateDerivatives<Dimension>& derivs,
           typename Dimension::Scalar /*currentTime*/) const {

  const double tiny = std::numeric_limits<double>::epsilon();

  // Get some useful fluid variables from the DataBase.
  const auto  mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  const auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity(),
                                                         this->requireOverlapConnectivity(),
                                                         this->requireIntersectionConnectivity());
  const auto  pairs = connectivityMap.nodePairList();

  // Initialize the return value to some impossibly high value.
  auto minDt = make_pair(std::numeric_limits<double>::max(), string());

  // Define a function for computing the velocity divergence (different for curvilinear coordinates)
  const auto Fdiv = (GeometryRegistrar::coords() == CoordinateType::Spherical ? +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi[0] + 2.0*veli[0]*safeInv(posi[0]); } :
                     GeometryRegistrar::coords() == CoordinateType::RZ        ? +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi.Trace() + veli[1]*safeInv(posi[1]); } :
                                                                                +[](const Tensor& DvDxi, const Vector& posi, const Vector& veli) { return DvDxi.Trace(); });

  // Loop over every fluid node.
  for (auto [nodeListi, fluidNodeListPtr]: enumerate(dataBase.fluidNodeListPtrs())) {
    const auto nPerh = fluidNodeListPtr->nodesPerSmoothingScale();
    CHECK(nPerh > 0.0);

    // Walk all the nodes in this FluidNodeList.
    const auto ni = connectivityMap.numNodes(nodeListi);
    const auto rank = Process::getRank();
#pragma omp parallel
    {
      auto minDt_local = minDt;
      auto DTrank_local = mDTrank;
      auto DTNodeList_local = mDTNodeList;
      auto DTnode_local = mDTnode;
      auto DTreason_local = mDTreason;
#pragma omp for
      for (auto k = 0u; k < ni; ++k) {
        const auto i = connectivityMap.ithNode(nodeListi, k);

        // If this node is masked, don't worry about it.
        if (mask(nodeListi, i) == 1) {

          // Get this nodes minimum characteristic smoothing scale.
          
          const auto& Hi = H(nodeListi, i);
          const auto& vi = velocity(nodeListi, i);
          CHECK2(Hi.Determinant() >  0.0,
                 "Bad H tensor : " << Hi << " : " << fluidNodeListPtr->name() << " " << i << " " << fluidNodeListPtr->firstGhostNode());
          const Scalar nodeScalei = 1.0/Hi.eigenValues().maxElement()/nPerh;

          // Velocity divergence limit.
          const auto divVelocity = Fdiv(DvDx(nodeListi, i), position(nodeListi, i), velocity(nodeListi, i));
          const auto divvDt = 1.0/(std::abs(divVelocity) + tiny);
          if (divvDt < minDt_local.first) {
            minDt_local = make_pair(divvDt, ("Velocity divergence limit: dt = " + to_string(divvDt) + "\n" +
                                             "                 div velocity = " + to_string(divVelocity) + "\n" +
                                             "                     material = " + fluidNodeListPtr->name() + "\n" +
                                             "        (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                             "                   @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "velocity divergence";
          }

          // Total acceleration limit.
          const auto vmagi = vi.magnitude();
          const auto dtAcc = 0.1*std::max(nodeScalei/(vmagi + tiny), vmagi/(DvDt(nodeListi, i).magnitude() + tiny));
          if (dtAcc < minDt_local.first) {
            minDt_local = make_pair(dtAcc, ("Total acceleration limit: dt = " + to_string(dtAcc) + "\n" + 
                                            "              |acceleration| = " + to_string(DvDt(nodeListi, i).magnitude()) + "\n" +
                                            "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                    material = " + fluidNodeListPtr->name() + "\n" +
                                            "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                            "                  @ position = " + vec_to_string(position(nodeListi, i))));
            DTrank_local = rank;
            DTNodeList_local = nodeListi;
            DTnode_local = i;
            DTreason_local = "acceleration";
          }

          // If requested, limit against the absolute velocity.
          if (useVelocityMagnitudeForDt()) {
            const auto velDt = nodeScalei/(velocity(nodeListi, i).magnitude() + 1.0e-10);
            if (velDt < minDt_local.first) {
              minDt_local = make_pair(velDt, ("Velocity magnitude limit: dt = " + to_string(velDt) + "\n" +
                                              "                        |vi| = " + to_string(velocity(nodeListi, i).magnitude()) + "\n" +
                                              "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                              "                    material = " + fluidNodeListPtr->name() + "\n" +
                                              "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                              "                  @ position = " + vec_to_string(position(nodeListi, i))));
              DTrank_local = rank;
              DTNodeList_local = nodeListi;
              DTnode_local = i;
              DTreason_local = "velocity magnitude";
            }
          }
        }
      }

      // Maximum pair-wise velocity difference limit
      const auto npairs = pairs.size();
#pragma omp for
      for (auto k = 0u; k < npairs; ++k) {
        const auto i = pairs[k].i_node;
        const auto j = pairs[k].j_node;
        const auto nodeListi = pairs[k].i_list;
        const auto nodeListj = pairs[k].j_list;

        const auto& vi = velocity(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  nodeScalei = 1.0/Hi.eigenValues().maxElement()/nPerh;

        const auto& vj = velocity(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  nodeScalej = 1.0/Hj.eigenValues().maxElement()/nPerh;

        const auto  vij = vi - vj;
        const auto  dtVelDiff = std::min(nodeScalei, nodeScalej)*safeInvVar(vij.magnitude(), tiny);
        if (dtVelDiff < minDt_local.first) {
          minDt_local = make_pair(dtVelDiff, ("Pairwise velocity difference limit: dt = " + to_string(dtVelDiff) + "\n" + 
                                              "                              material = " + fluidNodeListPtr->name() + "\n" +
                                              "                  (nodeListi, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(rank) + ")\n" +
                                              "                  (nodeListj, i, rank) = (" + to_string(nodeListj) + " " + to_string(j) + " " + to_string(rank) + ")\n" +
                                              "                   @ pos(nodeListi, i) = " + vec_to_string(position(nodeListi, i)) + "\n" +
                                              "                   @ pos(nodeListj, j) = " + vec_to_string(position(nodeListj, j)) + "\n" +
                                              "                                   vij = " + to_string(vij.magnitude()) + "\n" +
                                              "                            nodeScalei = " + to_string(nodeScalei) + "\n" +
                                              "                            nodeScalej = " + to_string(nodeScalej)));
          DTrank_local = rank;
          DTNodeList_local = nodeListi;
          DTnode_local = i;
          DTreason_local = "pairwise velocity difference";
        }
      }

#pragma omp critical
      {
        if (minDt_local.first < minDt.first) {
          minDt = minDt_local;
          mDTrank = DTrank_local;
          mDTNodeList = DTNodeList_local;
          mDTnode = DTnode_local;
          mDTreason = DTreason_local;
        }
      }
    }
  }

  // Scale by the cfl safety factor.
  minDt.first *= cfl();
  return minDt;
}

}
