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
#include "RK/computeVoronoiVolume.hh"
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
                   const TableKernel<Dimension>& W):
  SmoothingScaleBase<Dimension>(HUpdate),
  mWT(W),
  mZerothMoment(FieldStorageType::CopyFields),
  mFirstMoment(FieldStorageType::CopyFields),
  mSecondMoment(FieldStorageType::CopyFields),
  mCellSecondMoment(FieldStorageType::CopyFields),
  mCells(FieldStorageType::CopyFields),
  mDeltaCentroid(FieldStorageType::CopyFields) {
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
  dataBase.resizeFluidFieldList(mCells, FacetedVolume(), HydroFieldNames::cells, false);
  dataBase.resizeFluidFieldList(mDeltaCentroid, Vector::zero, "delta centroid", false);
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

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  CHECK(position.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(mass.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DHDt = derivs.fields(IncrementBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivs.fields(ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  massZerothMoment = derivs.fields(HydroFieldNames::massZerothMoment, 0.0);
  auto  massFirstMoment = derivs.fields(HydroFieldNames::massFirstMoment, Vector::zero);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(massZerothMoment.size() == numNodeLists);
  CHECK(massFirstMoment.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {
    // Thread private scratch variables
    bool sameMatij;
    int i, j, nodeListi, nodeListj;
    Scalar mi, mj, rhoi, rhoj, WSPHi, WSPHj, etaMagi, etaMagj, fweightij;
    Vector rij, etai, etaj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto massZerothMoment_thread = massZerothMoment.threadCopy(threadStack);
    auto massFirstMoment_thread = massFirstMoment.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      mi = mass(nodeListi, i);
      rhoi = massDensity(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);

      auto& massZerothMomenti = massZerothMoment_thread(nodeListi, i);
      auto& massFirstMomenti = massFirstMoment_thread(nodeListi, i);

      // Get the state for node j
      mj = mass(nodeListj, j);
      rhoj = massDensity(nodeListj, j);
      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);

      auto& massZerothMomentj = massZerothMoment_thread(nodeListj, j);
      auto& massFirstMomentj = massFirstMoment_thread(nodeListj, j);

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

      // Moments of the node distribution -- used for the ideal H calculation.
      fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
      massZerothMomenti +=     fweightij*WSPHi;
      massZerothMomentj += 1.0/fweightij*WSPHj;
      massFirstMomenti -=     fweightij*WSPHi*etai;
      massFirstMomentj += 1.0/fweightij*WSPHj*etaj;
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
  const auto& cm = dataBase.connectivityMap();
  const auto  pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  auto        H = state.fields(HydroFieldNames::H, SymTensor::zero);
  auto        Hideal = derivs.fields(ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);

  // Pair connectivity
  const auto& pairs = cm.nodePairList();
  const auto  npairs = pairs.size();

  // Compute the current Voronoi cells
  FieldList<Dimension, SymTensor> D;
  vector<Boundary<Dimension>*> boundaries(this->boundaryBegin(), this->boundaryEnd());
  auto vol = mass/rho;
  auto surfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  auto etaVoidPoints = dataBase.newFluidFieldList(vector<Vector>(), "etaVoidPoints");
  FieldList<Dimension, vector<CellFaceFlag>> cellFaceFlags;
  computeVoronoiVolume(pos, H, cm, D,
                       vector<FacetedVolume>(),          // facetedBoundaries
                       vector<vector<FacetedVolume>>(),  // holes
                       boundaries,
                       vol,                              // Use volume as weight
                       surfacePoint,
                       vol,
                       mDeltaCentroid,
                       etaVoidPoints,
                       mCells,
                       cellFaceFlags); 

  // Compute the second moments for the Voronoi cells
  const auto numNodeLists = dataBase.numFluidNodeLists();
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = mCells[k]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      mCellSecondMoment(k,i) = polySecondMoment(mCells(k,i), pos(k,i));
    }
  }

  // Apply boundary conditions to the cell second moment
  for (auto* boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->applyFieldListGhostBoundary(mCellSecondMoment);
  for (auto* boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();

  // Sum the net moments at each point
  mZerothMoment = 0.0;
  mSecondMoment = SymTensor::zero;
#pragma omp parallel
  {
    // Thread private scratch variables
    bool sameMatij;
    int i, j, nodeListi, nodeListj;
    Scalar mi, mj, rhoi, rhoj, WSPHi, WSPHj, etaMagi, etaMagj, fweightij;
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

      // Increment the moments for the pair
      fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
      massZerothMomenti +=     fweightij * WSPHi;
      massZerothMomentj += 1.0/fweightij * WSPHj;
      massSecondMomenti +=                 WSPHi*WSPHi * mCellSecondMoment(nodeListj, j);
      massSecondMomentj += 1.0/fweightij * WSPHj*WSPHj * mCellSecondMoment(nodeListi, i);
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
  }   // OpenMP parallel region

  // Apply boundary conditions to the moments
  for (auto* boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) {
    boundaryPtr->applyFieldListGhostBoundary(mZerothMoment);
    boundaryPtr->applyFieldListGhostBoundary(mSecondMoment);
  }
  for (auto* boundaryPtr: range(this->boundaryBegin(), this->boundaryEnd())) boundaryPtr->finalizeGhostBoundary();

  // Now we have the moments, so we can loop over the points and set our new H
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto& nodeList = mass[k]->nodeList();
    const auto  hminInv = safeInvVar(nodeList.hmin());
    const auto  hmaxInv = safeInvVar(nodeList.hmax());
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();
    const auto  n = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      auto&       Hi = H(k,i);
      auto&       Hideali = Hideal(k,i);
      auto        massZerothMomenti = mZerothMoment(k,i);
      const auto& massSecondMomenti = mSecondMoment(k,i);
  
      // Complete the zeroth moment
      massZerothMomenti = Dimension::rootnu(max(0.0, massZerothMomenti));

      // Find the new normalized target shape
      auto T = massSecondMomenti.sqrt();
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
    }
  }
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
