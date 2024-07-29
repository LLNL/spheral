//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScale
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 15:01:13 PDT 2005
//----------------------------------------------------------------------------//
#include "SmoothingScale/ASPHSmoothingScale.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/range.hh"
#include "Utilities/Timer.hh"

#include <cmath>
#include <vector>

namespace Spheral {

using std::min;
using std::max;
using std::abs;
using std::vector;

namespace {

//------------------------------------------------------------------------------
// DH/Dt per dimension
//------------------------------------------------------------------------------
// 1-D case same as SPH.
inline
Dim<1>::SymTensor
smoothingScaleDerivative(const Dim<1>::SymTensor& H,
                         const Dim<1>::Tensor& DvDx) {
  return -H*DvDx.Trace();
}

// 2-D ASPH tensor evolution.
inline
Dim<2>::SymTensor
smoothingScaleDerivative(const Dim<2>::SymTensor& H,
                         const Dim<2>::Tensor& DvDx) {
  REQUIRE(H.Trace() > 0.0);
  const auto thetaDot = (H.xx()*DvDx.xy() - H.yy()*DvDx.yx() - H.yx()*(DvDx.xx() - DvDx.yy()))/H.Trace();
  Dim<2>::SymTensor result;
  result.xx(H.yx()*(thetaDot - DvDx.yx()) - H.xx()*DvDx.xx());
  result.xy(-(H.xx()*thetaDot + H.yx()*DvDx.xx() + H.yy()*DvDx.yx()));
  result.yy(-H.yx()*(thetaDot + DvDx.xy()) - H.yy()*DvDx.yy());
  return result;
}

// 3-D ASPH tensor evolution.
inline
Dim<3>::SymTensor
smoothingScaleDerivative(const Dim<3>::SymTensor& H,
                         const Dim<3>::Tensor& DvDx) {
  REQUIRE(H.Trace() > 0.0);
  const auto AA = H.xx()*DvDx.xy() - H.xy()*(DvDx.xx() - DvDx.yy()) + H.xz()*DvDx.zy() - H.yy()*DvDx.yx() - H.yz()*DvDx.zx();
  const auto BB = H.xx()*DvDx.xz() + H.xy()*DvDx.yz() - H.xz()*(DvDx.xx() - DvDx.zz()) - H.yz()*DvDx.yx() - H.zz()*DvDx.zx();
  const auto CC = H.xy()*DvDx.xz() + H.yy()*DvDx.yz() - H.yz()*(DvDx.yy() - DvDx.zz()) - H.xz()*DvDx.xy() - H.zz()*DvDx.zy();
  const auto thpt = H.yy() + H.zz();
  const auto Ga = (H.xx() + H.yy())*thpt - H.xz()*H.xz();
  const auto Gb = (H.yy() + H.zz())*H.yz() + H.xy()*H.xz();
  const auto Gc = (H.xx() + H.zz())*thpt - H.xy()*H.xy();
  const auto Gd = thpt*AA + H.xz()*CC;
  const auto Ge = thpt*BB - H.xy()*CC;
  const auto ack = 1.0/(Ga*Gc - Gb*Gb);
  const auto Gdot = (Gc*Gd - Gb*Ge)*ack;
  const auto Tdot = (Gb*Gd - Ga*Ge)*ack;
  const auto Phidot = (H.xz()*Gdot + H.xy()*Tdot + CC)/thpt;
  Dim<3>::SymTensor result;
  result.xx(-H.xx()*DvDx.xx() + H.xy()*(Gdot - DvDx.yx()) - H.xz()*(Tdot + DvDx.zx()));
  result.xy(H.yy()*Gdot - H.yz()*Tdot - H.xx()*DvDx.xy() - H.xy()*DvDx.yy() - H.xz()*DvDx.zy());
  result.xz(H.yz()*Gdot - H.zz()*Tdot - H.xx()*DvDx.xz() - H.xy()*DvDx.yz() - H.xz()*DvDx.zz());
  result.yy(H.yz()*(Phidot - DvDx.zy()) - H.xy()*(Gdot + DvDx.xy()) - H.yy()*DvDx.yy());
  result.yz(H.xy()*Tdot - H.yy()*Phidot - H.xz()*DvDx.xy() - H.yz()*DvDx.yy() - H.zz()*DvDx.zy());
  result.zz(H.xz()*(Tdot - DvDx.xz()) - H.yz()*(Phidot + DvDx.yz()) - H.zz()*DvDx.zz());
  return result;
}
  
//------------------------------------------------------------------------------
// Compute the second moment about the give position for a polytope
//------------------------------------------------------------------------------
// 1D -- nothing to do
inline
Dim<1>::SymTensor
polySecondMoment(const Dim<1>::FacetedVolume& poly,
                 const Dim<1>::Vector& center) {
  return Dim<1>::SymTensor(1);
}

// 2D -- we can use the knowledge that the vertices in a 
inline
Dim<2>::SymTensor
polySecondMoment(const Dim<2>::FacetedVolume& poly,
                 const Dim<2>::Vector& center) {
  Dim<2>::SymTensor result;
  const auto& facets = poly.facets();
  for (const auto& f: facets) {
    const auto v1 = f.point1() - center;
    const auto v2 = f.point2() - center;
    const auto thpt = std::abs(v1.x()*v2.y() - v2.x()*v1.y())/12.0;
    result[0] += (v1.x()*v1.x() + v1.x()*v2.x() + v2.x()*v2.x())*thpt;
    result[1] += (v1.x()*v1.y() + v2.x()*v2.y() + 0.5*(v2.x()*v1.y() + v1.x()*v2.y()))*thpt;
    result[2] += (v1.y()*v1.y() + v1.y()*v2.y() + v2.y()*v2.y())*thpt;
  }
  return result;
}

inline
Dim<3>::SymTensor
polySecondMoment(const Dim<3>::FacetedVolume& poly,
                 const Dim<3>::Vector& center) {
  VERIFY2(false, "Implement me!");
  return Dim<3>::SymTensor();
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScale<Dimension>::
ASPHSmoothingScale(const HEvolutionType HUpdate,
                   const TableKernel<Dimension>& W,
                   const Scalar fHourGlass):
  SmoothingScaleBase<Dimension>(HUpdate),
  mfHourGlass(fHourGlass),
  mWT(W),
  mZerothMoment(FieldStorageType::CopyFields),
  mFirstMoment(FieldStorageType::CopyFields),
  mSecondMoment(FieldStorageType::CopyFields),
  mCellSecondMoment(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  // Make sure our FieldLists are correctly sized.
  SmoothingScaleBase<Dimension>::initializeProblemStartup(dataBase);
  dataBase.resizeFluidFieldList(mZerothMoment, 0.0, HydroFieldNames::massZerothMoment, false);
  dataBase.resizeFluidFieldList(mFirstMoment, Vector::zero, HydroFieldNames::massFirstMoment, false);
  dataBase.resizeFluidFieldList(mSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mCellSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment + " cells", false);
}

//------------------------------------------------------------------------------
// On problem start up (following above), we need initialize the cell geometries
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // // Grab our state
  // const auto numNodeLists = dataBase.numFluidNodeLists();
  // const auto  pos = state.fields(HydroFieldNames::position, Vector::zero);
  // const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  // const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  // const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);

  // // Connectivity
  // dataBase.updateConnectivityMap(false, false, false);
  // const auto& cm = dataBase.connectivityMap();

  // // Compute the current Voronoi cells
  // FieldList<Dimension, SymTensor> D;
  // vector<Boundary<Dimension>*> boundaries(this->boundaryBegin(), this->boundaryEnd());
  // auto vol = mass/rho;
  // auto surfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  // auto etaVoidPoints = dataBase.newFluidFieldList(vector<Vector>(), "etaVoidPoints");
  // FieldList<Dimension, vector<CellFaceFlag>> cellFaceFlags;
  // computeVoronoiVolume(pos, H, cm, D,
  //                      vector<FacetedVolume>(),          // facetedBoundaries
  //                      vector<vector<FacetedVolume>>(),  // holes
  //                      boundaries,
  //                      FieldList<Dimension, Scalar>(),   // weight
  //                      surfacePoint,
  //                      vol,
  //                      mDeltaCentroid,
  //                      etaVoidPoints,
  //                      mCells,
  //                      cellFaceFlags); 
}

//------------------------------------------------------------------------------
// Register state
// Override the normal SmoothingScaleBase version since we only do the idealH
// update in the finalize step at the end of advancement step (due to expense).
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  const auto Hupdate = this->HEvolution();
  auto Hfields = dataBase.fluidHfield();
  const auto numFields = Hfields.numFields();
  for (auto k = 0u; k < numFields; ++k) {
    auto& Hfield = *Hfields[k];
    const auto& nodeList = Hfield.nodeList();
    const auto hmaxInv = 1.0/nodeList.hmax();
    const auto hminInv = 1.0/nodeList.hmin();
    switch (Hupdate) {
      case HEvolutionType::IntegrateH:
      case HEvolutionType::IdealH:
        state.enroll(Hfield, make_policy<IncrementBoundedState<Dimension, SymTensor, Scalar>>(hmaxInv, hminInv));
        break;

      case HEvolutionType::FixedH:
        state.enroll(Hfield);
        break;

       default:
         VERIFY2(false, "ASPHSmoothingScale ERROR: Unknown Hevolution option ");
    }
  }
}

//------------------------------------------------------------------------------
// Register derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  SmoothingScaleBase<Dimension>::registerDerivatives(dataBase, derivs);
  derivs.enroll(mZerothMoment);
  derivs.enroll(mFirstMoment);
  derivs.enroll(mSecondMoment);
}

//------------------------------------------------------------------------------
// Time derivative of the smoothing scale.
// We depend on a previous package evaluating the velocity gradient (DvDx)
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("ASPHSmoothingScaleDerivs");

  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto P = state.fields(HydroFieldNames::pressure, 0.0);
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  CHECK(position.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(mass.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DHDt = derivs.fields(IncrementBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivs.fields(ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  massZerothMoment = derivs.fields(HydroFieldNames::massZerothMoment, 0.0);
  auto  massFirstMoment = derivs.fields(HydroFieldNames::massFirstMoment, Vector::zero);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(massZerothMoment.size() == numNodeLists);
  CHECK(massFirstMoment.size() == numNodeLists);

  // // Check if we're using a compatible discretization for the momentum & energy
  // auto& pairAccelerations = derivs.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  // const bool compatibleEnergy = (pairAccelerations.size() == npairs);
  // const bool useHourGlass = (mCells.size() == numNodeLists and mfHourGlass > 0.0);

#pragma omp parallel
  {
    // Thread private scratch variables
    bool sameMatij;
    int i, j, nodeListi, nodeListj;
    Scalar mi, mj, rhoi, rhoj, Pi, Pj, Pij, WSPHi, WSPHj, etaMagi, etaMagj, fweightij;
    Scalar Wi, Wj;
    Vector rij, etai, etaj, gradWi, gradWj;
    SymTensor psiij;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto massZerothMoment_thread = massZerothMoment.threadCopy(threadStack);
    auto massFirstMoment_thread = massFirstMoment.threadCopy(threadStack);
    auto DvDt_thread = DvDt.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      mi = mass(nodeListi, i);
      rhoi = massDensity(nodeListi, i);
      // Pi = P(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);

      auto& massZerothMomenti = massZerothMoment_thread(nodeListi, i);
      auto& massFirstMomenti = massFirstMoment_thread(nodeListi, i);
      // auto& DvDti = DvDt_thread(nodeListi, i);

      // Get the state for node j
      mj = mass(nodeListj, j);
      rhoj = massDensity(nodeListj, j);
      // Pj = P(nodeListj, j);
      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);

      auto& massZerothMomentj = massZerothMoment_thread(nodeListj, j);
      auto& massFirstMomentj = massFirstMoment_thread(nodeListj, j);
      // auto& DvDtj = DvDt_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      sameMatij = (nodeListi == nodeListj); // and fragIDi == fragIDj);

      // Node displacement.
      rij = ri - rj;
      etai = Hi*rij;
      etaj = Hj*rij;
      etaMagi = etai.magnitude();
      etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      WSPHi = mWT.kernelValueSPH(etaMagi);
      WSPHj = mWT.kernelValueSPH(etaMagj);
      gradWi = mWT.gradValue(etaMagi, Hi.Determinant()) * Hi*etai*safeInvVar(etaMagi);
      gradWj = mWT.gradValue(etaMagj, Hj.Determinant()) * Hj*etaj*safeInvVar(etaMagj);

      // Moments of the node distribution -- used for the ideal H calculation.
      fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
      psiij = rij.unitVector().selfdyad();
      massZerothMomenti +=     fweightij*WSPHi;
      massZerothMomentj += 1.0/fweightij*WSPHj;
      massFirstMomenti -=     fweightij*WSPHi*etai;
      massFirstMomentj += 1.0/fweightij*WSPHj*etaj;

      // // Add term to fight pairing instability with high-aspect ratio points
      // if (useHourGlass) {
      //   const auto centi = mDeltaCentroid(nodeListi, i); // mCells(nodeListi, i).centroid();
      //   const auto centj = mDeltaCentroid(nodeListj, j); // mCells(nodeListj, j).centroid();
      //   const auto  cij = centi - centj;
      //   const auto  cijMag = cij.magnitude();
      //   CHECK(cijMag > 0.0);
      //   const auto  chat = cij/cijMag;
      //   Pij = mfHourGlass * max(abs(Pi), abs(Pj)) * (1.0 - min(1.0, abs(rij.dot(chat))/cijMag));
      //   CHECK(Pij >= 0.0);
      //   const auto deltaDvDt = Pij/(rhoi*rhoi)*gradWi + Pij/(rhoj*rhoj)*gradWj;
      //   DvDti -= mj*deltaDvDt;
      //   DvDtj += mi*deltaDvDt;
      //   if (compatibleEnergy) pairAccelerations[kk] -= mj*deltaDvDt;
      // }
    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
  }   // OpenMP parallel region

  // Finish up the derivatives now that we've walked all pairs
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hminInv = safeInvVar(nodeList.hmin());
    const auto  hmaxInv = safeInvVar(nodeList.hmax());
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& Hi = H(nodeListi, i);
      const auto& DvDxi = DvDx(nodeListi, i);

      auto& massZerothMomenti = massZerothMoment(nodeListi, i);
      // const auto& massFirstMomenti = massFirstMoment(nodeListi, i);
      // const auto& massSecondMomenti = massSecondMoment(nodeListi, i);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      massZerothMomenti = Dimension::rootnu(max(0.0, massZerothMomenti));

      // Time derivative of H
      DHDt(nodeListi, i) = smoothingScaleDerivative(Hi, DvDxi);

      // Determine the current effective number of nodes per smoothing scale.
      const auto currentNodesPerSmoothingScale = (fuzzyEqual(massZerothMomenti, 0.0) ?  // Is this node isolated (no neighbors)?
                                                  0.5*nPerh :
                                                  mWT.equivalentNodesPerSmoothingScale(massZerothMomenti));
      CHECK2(currentNodesPerSmoothingScale > 0.0, "Bad estimate for nPerh effective from kernel: " << currentNodesPerSmoothingScale);

      // The ratio of the desired to current nodes per smoothing scale.
      const auto s = std::min(4.0, std::max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30)));
      CHECK(s > 0.0);

      // Now determine how to scale the current H to the desired value.
      // We only scale H at this point, not try to change the shape.
      const auto a = (s < 1.0 ? 
                      0.4*(1.0 + s*s) :
                      0.4*(1.0 + 1.0/(s*s*s)));
      CHECK(1.0 - a + a*s > 0.0);
      Hideal(nodeListi, i) = std::max(hmaxInv, std::min(hminInv, Hi / (1.0 - a + a*s)));
    }
  }
  TIME_END("ASPHSmoothingScaleDerivs");
}

//------------------------------------------------------------------------------
// Finalize at the end of the step.
// This is where we compute the Voronoi cell geometry and use it to set our
// second moments and new H shape.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
finalize(const Scalar time, 
         const Scalar dt,
         DataBase<Dimension>& dataBase, 
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Grab our state
  const auto numNodeLists = dataBase.numFluidNodeLists();
  const auto& cm = dataBase.connectivityMap();
  auto        pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  vel = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  cells = state.fields(HydroFieldNames::cells, FacetedVolume());
  const auto  surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto        H = state.fields(HydroFieldNames::H, SymTensor::zero);
  auto        Hideal = derivs.fields(ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);

  // Pair connectivity
  const auto& pairs = cm.nodePairList();
  const auto  npairs = pairs.size();

  // Compute the second moments for the Voronoi cells
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = cells[k]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      mCellSecondMoment(k,i) = polySecondMoment(cells(k,i), pos(k,i)).sqrt();
    }
  }

  // Apply boundary conditions to the cell second moments
  for (auto* boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->applyFieldListGhostBoundary(mCellSecondMoment);
    boundaryPtr->finalizeGhostBoundary();
  }

// //   // Prepare RK correction terms
// //   FieldList<Dimension, Scalar> m0 = dataBase.newFluidFieldList(0.0, "m0");
// //   FieldList<Dimension, Vector> m1 = dataBase.newFluidFieldList(Vector::zero, "m1");
// //   FieldList<Dimension, SymTensor> m2 = dataBase.newFluidFieldList(SymTensor::zero, "m2");
// //   FieldList<Dimension, Scalar> A = dataBase.newFluidFieldList(0.0, "A");
// //   FieldList<Dimension, Vector> B = dataBase.newFluidFieldList(Vector::zero, "B");
// // #pragma omp parallel
// //   {
// //     // Thread private scratch variables
// //     bool sameMatij;
// //     int i, j, nodeListi, nodeListj;
// //     Scalar mi, mj, rhoi, rhoj, WSPHi, WSPHj, etaMagi, etaMagj, fweightij;
// //     Vector rij, etai, etaj;

// //     typename SpheralThreads<Dimension>::FieldListStack threadStack;
// //     auto m0_thread = m0.threadCopy(threadStack);
// //     auto m1_thread = m1.threadCopy(threadStack);
// //     auto m2_thread = m2.threadCopy(threadStack);

// // #pragma omp for
// //     for (auto kk = 0u; kk < npairs; ++kk) {
// //       i = pairs[kk].i_node;
// //       j = pairs[kk].j_node;
// //       nodeListi = pairs[kk].i_list;
// //       nodeListj = pairs[kk].j_list;

// //       // State for node i
// //       mi = mass(nodeListi, i);
// //       rhoi = rho(nodeListi, i);
// //       const auto& ri = pos(nodeListi, i);
// //       const auto& Hi = H(nodeListi, i);
// //       auto& m0i = m0_thread(nodeListi, i);
// //       auto& m1i = m1_thread(nodeListi, i);
// //       auto& m2i = m2_thread(nodeListi, i);

// //       // Get the state for node j
// //       mj = mass(nodeListj, j);
// //       rhoj = rho(nodeListj, j);
// //       const auto& rj = pos(nodeListj, j);
// //       const auto& Hj = H(nodeListj, j);
// //       auto& m0j = m0_thread(nodeListj, j);
// //       auto& m1j = m1_thread(nodeListj, j);
// //       auto& m2j = m2_thread(nodeListj, j);

// //       // Flag if this is a contiguous material pair or not.
// //       sameMatij = (nodeListi == nodeListj); // and fragIDi == fragIDj);

// //       // Node displacement.
// //       rij = ri - rj;
// //       etai = Hi*rij;
// //       etaj = Hj*rij;
// //       etaMagi = etai.magnitude();
// //       etaMagj = etaj.magnitude();
// //       CHECK(etaMagi >= 0.0);
// //       CHECK(etaMagj >= 0.0);

// //       // Symmetrized kernel weight and gradient.
// //       WSPHi = mWT.kernelValueSPH(etaMagi);
// //       WSPHj = mWT.kernelValueSPH(etaMagj);

// //       // Sum the moments
// //       fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
// //       m0i +=     fweightij * WSPHi;
// //       m0j += 1.0/fweightij * WSPHj;
// //       m1i +=     fweightij * WSPHi*rij;
// //       m1j -= 1.0/fweightij * WSPHj*rij;
// //       m2i +=     fweightij * WSPHi*rij.selfdyad();
// //       m2j += 1.0/fweightij * WSPHj*rij.selfdyad();
// //     }

// //     // Reduce the thread values to the master.
// //     threadReduceFieldLists<Dimension>(threadStack);
// //   }   // OpenMP parallel region
      
// //   // Compute the corrections
// //   for (auto k = 0u; k < numNodeLists; ++k) {
// //     const auto& nodeList = mass[k]->nodeList();
// //     const auto  n = nodeList.numInternalNodes();
// // #pragma omp parallel for
// //     for (auto i = 0u; i < n; ++i) {
// //       A(k,i) = 1.0/(m0(k,i) - m2(k,i).Inverse().dot(m1(k,i)).dot(m1(k,i)));
// //       B(k,i) = -m2(k,i).Inverse().dot(m1(k,i));
// //     }
// //   }

  // Sum the net moments at each point
  mZerothMoment = 0.0;
  mSecondMoment = SymTensor::zero;
#pragma omp parallel
  {
    // Thread private scratch variables
    bool sameMatij;
    int i, j, nodeListi, nodeListj;
    Scalar mi, mj, rhoi, rhoj, WSPHi, WSPHj, WRKi, WRKj, etaMagi, etaMagj, fweightij;
    Scalar Wi, Wj;
    Vector rij, etai, etaj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto massZerothMoment_thread = mZerothMoment.threadCopy(threadStack);
    auto massSecondMoment_thread = mSecondMoment.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // State for node i
      mi = mass(nodeListi, i);
      rhoi = rho(nodeListi, i);
      const auto& ri = pos(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      auto& massZerothMomenti = massZerothMoment_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);

      // Get the state for node j
      mj = mass(nodeListj, j);
      rhoj = rho(nodeListj, j);
      const auto& rj = pos(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      auto& massZerothMomentj = massZerothMoment_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      sameMatij = (nodeListi == nodeListj); // and fragIDi == fragIDj);

      // Node displacement.
      rij = ri - rj;
      etai = Hi*rij;
      etaj = Hj*rij;
      etaMagi = etai.magnitude();
      etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      WSPHi = mWT.kernelValueSPH(etaMagi);
      WSPHj = mWT.kernelValueSPH(etaMagj);
      // Wi = mWT.kernelValue(etaMagi, 1.0);
      // Wj = mWT.kernelValue(etaMagj, 1.0);
      // WRKi = WSPHi * A(nodeListi, i)*(1.0 - B(nodeListi, i).dot(rij));
      // WRKj = WSPHj * A(nodeListj, j)*(1.0 + B(nodeListj, j).dot(rij));

      // Increment the moments for the pair
      fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
      massZerothMomenti +=     fweightij * WSPHi;
      massZerothMomentj += 1.0/fweightij * WSPHj;
      massSecondMomenti +=                 WSPHi * mCellSecondMoment(nodeListj, j);
      massSecondMomentj += 1.0/fweightij * WSPHj * mCellSecondMoment(nodeListi, i);
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
  }   // OpenMP parallel region

  // // Apply boundary conditions to the moments
  // for (auto* boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
  //   boundaryPtr->applyFieldListGhostBoundary(mZerothMoment);
  //   boundaryPtr->applyFieldListGhostBoundary(mSecondMoment);
  // }
  // for (auto* boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();

  // Now we have the moments, so we can loop over the points and set our new H
  const auto W0 = mWT.kernelValue(0.0, 1.0);
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto& nodeList = mass[k]->nodeList();
    const auto  hminInv = safeInvVar(nodeList.hmin());
    const auto  hmaxInv = safeInvVar(nodeList.hmax());
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();
    const auto  n = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      auto&       ri = pos(k,i);
      auto&       Hi = H(k,i);
      auto&       Hideali = Hideal(k,i);
      auto        massZerothMomenti = mZerothMoment(k,i);
      auto&       massSecondMomenti = mSecondMoment(k,i);
  
      // Complete the zeroth moment
      massZerothMomenti = Dimension::rootnu(max(0.0, massZerothMomenti));

      // // Complete the second moment
      // massSecondMomenti += W0 * polySecondMoment(mCells(k,i), ri).sqrt();

      // Find the new normalized target shape
      auto T = massSecondMomenti; // .sqrt();
      {
        const auto detT = T.Determinant();
        if (fuzzyEqual(detT, 0.0)) {
          T = SymTensor::one;
        } else {
          T /= Dimension::rootnu(detT);
        }
      }
      CHECK(fuzzyEqual(T.Determinant(), 1.0));
      T /= Dimension::rootnu(Hi.Determinant());   // T in units of length, now with same volume as the old Hinverse
      CHECK(fuzzyEqual(T.Determinant(), 1.0/Hi.Determinant()));
      
      // Determine the current effective number of nodes per smoothing scale.
      const auto currentNodesPerSmoothingScale = (fuzzyEqual(massZerothMomenti, 0.0) ?  // Is this node isolated (no neighbors)?
                                                  0.5*nPerh :
                                                  mWT.equivalentNodesPerSmoothingScale(massZerothMomenti));
      CHECK2(currentNodesPerSmoothingScale > 0.0, "Bad estimate for nPerh effective from kernel: " << currentNodesPerSmoothingScale);

      // The ratio of the desired to current nodes per smoothing scale.
      const auto s = std::min(4.0, std::max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30)));
      CHECK(s > 0.0);

      // // Determine the desired H determinant using our usual target nperh logic
      // auto fscale = 1.0;
      // for (auto j = 0u; j < Dimension::nDim; ++j) {
      //   eigenT.eigenValues[j] = std::max(eigenT.eigenValues[j], hminratio*Tmax);
      //   fscale *= eigenT.eigenValues[j];
      // }
      // CHECK(fscale > 0.0);
      // fscale = 1.0/Dimension::rootnu(fscale);

      // Now apply the desired volume scaling from the zeroth moment to fscale
      const auto a = (s < 1.0 ? 
                      0.4*(1.0 + s*s) :
                      0.4*(1.0 + 1.0/(s*s*s)));
      CHECK(1.0 - a + a*s > 0.0);
      T *= std::min(4.0, std::max(0.25, 1.0 - a + a*s));

      // Build the new H tensor
      // Hi = constructSymTensorWithBoundedDiagonal(fscale*eigenT.eigenValues, hmaxInv, hminInv);
      // Hi.rotationalTransform(eigenT.eigenVectors);
      Hi = T.Inverse();
      Hideali = Hi;                                 // To be consistent with SPH package behaviour

      // // If requested, move toward the cell centroid
      // if (mfHourGlass > 0.0 and surfacePoint(k,i) == 0) {
      //   const auto& vi = vel(k,i);
      //   const auto  ci = cs(k,i);
      //   const auto  vhat = vi*safeInv(vi.magnitude());   // goes to zero when velocity zero
      //   const auto  centi = cells(k,i).centroid();
      //   auto        dr = mfHourGlass*(centi - ri);
      //   dr = dr.dot(vhat) * vhat;
      //   // const auto  drmax = mfHourGlass*dt*vi.magnitude();
      //   const auto  drmax = mfHourGlass*dt*ci;
      //   // const auto  drmax = 0.5*dt*min(ci, vi.magnitude());
      //   const auto  drmag = dr.magnitude();
      //   dr *= min(1.0, drmax*safeInv(drmag));
      //   ri += dr;
      // }
    }
  }
}

//------------------------------------------------------------------------------
// Apply boundary conditions to the physics specific fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  SmoothingScaleBase<Dimension>::dumpState(file, pathName);
  file.write(mZerothMoment, pathName + "/zerothMoment");
  file.write(mFirstMoment, pathName + "/firstMoment");
  file.write(mSecondMoment, pathName + "/secondMoment");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  SmoothingScaleBase<Dimension>::restoreState(file, pathName);
  file.read(mZerothMoment, pathName + "/zerothMoment");
  file.read(mFirstMoment, pathName + "/firstMoment");
  file.read(mSecondMoment, pathName + "/secondMoment");
}

}
