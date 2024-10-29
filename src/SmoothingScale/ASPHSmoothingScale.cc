//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScale
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 15:01:13 PDT 2005
//----------------------------------------------------------------------------//
#include "SmoothingScale/ASPHSmoothingScale.hh"
#include "SmoothingScale/polySecondMoment.hh"
#include "SmoothingScale/IncrementASPHHtensor.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/range.hh"
#include "Utilities/Timer.hh"

#include <cmath>
#include <vector>

namespace Spheral {

using std::min;
using std::max;
using std::abs;
using std::vector;
using FastMath::pow2;

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
// Radial evolution/alignment specialized evolution
//------------------------------------------------------------------------------
// 1-D
inline
Dim<1>::SymTensor
radialEvolution(const Dim<1>::SymTensor& Hi,
                const Dim<1>::Vector& nhat,
                const Dim<1>::Scalar s,
                const Dim<1>::Scalar r0,
                const Dim<1>::Scalar r1) {
  return Hi / s;
}

// 2-D
inline
Dim<2>::SymTensor
radialEvolution(const Dim<2>::SymTensor& Hi,
                const Dim<2>::Vector& nhat,
                const Dim<2>::Scalar s,
                const Dim<2>::Scalar r0,
                const Dim<2>::Scalar r1) {
  const auto T = rotationMatrix(nhat).Transpose();
  const auto hev0 = Hi.eigenVectors();
  Dim<2>::SymTensor result;
  if (abs(hev0.eigenVectors.getColumn(0).dot(nhat)) > abs(hev0.eigenVectors.getColumn(1).dot(nhat))) {
    result(0,0) = hev0.eigenValues(0);
    result(1,1) = hev0.eigenValues(1);
  } else {
    result(0,0) = hev0.eigenValues(1);
    result(1,1) = hev0.eigenValues(0);
  }
  const auto fr = r1*safeInvVar(r0);
  CHECK(fr > 0.0);
  result(0,0) /= s;
  result(1,1) /= fr;
  result.rotationalTransform(T);
  return result;
}

// 3-D
inline
Dim<3>::SymTensor
radialEvolution(const Dim<3>::SymTensor& Hi,
                const Dim<3>::Vector& nhat,
                const Dim<3>::Scalar s,
                const Dim<3>::Scalar r0,
                const Dim<3>::Scalar r1) {
  const auto Tprinciple = rotationMatrix(nhat);
  const auto Tlab = Tprinciple.Transpose();
  auto result = Hi;
  result.rotationalTransform(Tprinciple);
  const auto fr = r1*safeInvVar(r0);
  CHECK(fr > 0.0);
  result(0,0) /= s;
  result(1,1) /= fr;
  result(1,2) /= fr;
  result(2,1) /= fr;
  result(2,2) /= fr;
  result.rotationalTransform(Tlab);
  return result;
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHSmoothingScale<Dimension>::
ASPHSmoothingScale(const HEvolutionType HUpdate,
                   const TableKernel<Dimension>& W,
                   const bool fixShape,
                   const bool radialOnly):
  SmoothingScaleBase<Dimension>(HUpdate),
  mWT(W),
  mZerothMoment(FieldStorageType::CopyFields),
  mSecondMoment(FieldStorageType::CopyFields),
  mCellSecondMoment(FieldStorageType::CopyFields),
  mRadius0(FieldStorageType::CopyFields),
  mHidealFilterPtr(std::make_shared<ASPHSmoothingScaleUserFilter<Dimension>>()),
  mRadialFunctorPtr(std::make_shared<ASPHRadialFunctor<Dimension>>()),
  mFixShape(fixShape),
  mRadialOnly(radialOnly) {
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
  dataBase.resizeFluidFieldList(mSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mCellSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment + " cells", false);
  if (mRadialOnly) dataBase.resizeFluidFieldList(mRadius0, 0.0, "Start of step radius", false);
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
  switch (Hupdate) {
  case HEvolutionType::IntegrateH:
  case HEvolutionType::IdealH:
    state.enroll(Hfields, make_policy<IncrementASPHHtensor<Dimension>>(mFixShape, mRadialOnly, mRadialFunctorPtr));
    break;

  case HEvolutionType::FixedH:
    state.enroll(Hfields);
    break;

  default:
    VERIFY2(false, "ASPHSmoothingScale ERROR: Unknown Hevolution option ");
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
  derivs.enroll(mSecondMoment);
  derivs.enroll(mCellSecondMoment);
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
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  CHECK(H.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);

  // Derivative FieldLists.
  auto DHDt = derivs.fields(IncrementBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  CHECK(DHDt.size() == numNodeLists);

  // Set the H time derivatives
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto& nodeList = H[k]->nodeList();
    const auto  ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto& Hi = H(k,i);
      const auto& DvDxi = DvDx(k,i);
      DHDt(k,i) = smoothingScaleDerivative(Hi, DvDxi);
    }
  }
  TIME_END("ASPHSmoothingScaleDerivs");
}

//------------------------------------------------------------------------------
// Initialize at the beginning of a timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHSmoothingScale<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  // If we're using the radial H scaling, take a snapshot of the initial radius of
  // each point.
  if (mRadialOnly) {
    const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
    const auto numFields = pos.numFields();
    for (auto k = 0u; k < numFields; ++k) {
      const auto n = pos[k]->numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        mRadius0(k,i) = mRadialFunctorPtr->radialCoordinate(k, i, pos(k,i));
      }
    }
  }
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

  // Notify any user filter object things are about to start
  mHidealFilterPtr->startFinalize(time, dt, dataBase, state, derivs);

  // If we're not using the IdealH algorithm we can save a lot of time...
  const auto Hupdate = this->HEvolution();
  if (Hupdate == HEvolutionType::IdealH) {
    CHECK(not (mFixShape and mRadialOnly));                // Can't do both simultaneously
    const auto voronoi = not (mFixShape or mRadialOnly);   // In these special cases we don't need the Voronoi second moment

    // Grab our state
    const auto numNodeLists = dataBase.numFluidNodeLists();
    const auto& cm = dataBase.connectivityMap();
    auto        pos = state.fields(HydroFieldNames::position, Vector::zero);
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
    const auto  cells = state.fields(HydroFieldNames::cells, FacetedVolume());
    const auto  surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
    auto        H = state.fields(HydroFieldNames::H, SymTensor::zero);
    auto        Hideal = derivs.fields(ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
    CHECK(pos.size() == numNodeLists);
    CHECK(mass.size() == numNodeLists);
    CHECK(rho.size() == numNodeLists);
    CHECK2((cells.size() == numNodeLists) or not voronoi, cells.size() << " " << voronoi << " " << mFixShape << " " << mRadialOnly);
    CHECK2((surfacePoint.size() == numNodeLists) or not voronoi, cells.size() << " " << voronoi << " " << mFixShape << " " << mRadialOnly);
    CHECK(H.size() == numNodeLists);
    CHECK(Hideal.size() == numNodeLists);

    // Pair connectivity
    const auto& pairs = cm.nodePairList();
    const auto  npairs = pairs.size();

    // Compute the second moments for the Voronoi cells
    if (voronoi) {
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
        // Wi = mWT.kernelValue(etaMagi, 1.0);
        // Wj = mWT.kernelValue(etaMagj, 1.0);
        // WRKi = WSPHi * A(nodeListi, i)*(1.0 - B(nodeListi, i).dot(rij));
        // WRKj = WSPHj * A(nodeListj, j)*(1.0 + B(nodeListj, j).dot(rij));

        // Increment the moments for the pair
        fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
        massZerothMomenti +=     fweightij * WSPHi;
        massZerothMomentj += 1.0/fweightij * WSPHj;
        if (voronoi) {
          if (surfacePoint(nodeListj, j) == 0) massSecondMomenti +=                 WSPHi * mCellSecondMoment(nodeListj, j);
          if (surfacePoint(nodeListi, i) == 0) massSecondMomentj += 1.0/fweightij * WSPHj * mCellSecondMoment(nodeListi, i);
        }
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
    // const auto W0 = mWT.kernelValue(0.0, 1.0);
    for (auto k = 0u; k < numNodeLists; ++k) {
      const auto& nodeList = mass[k]->nodeList();
      const auto  hminInv = safeInvVar(nodeList.hmin());
      const auto  hmaxInv = safeInvVar(nodeList.hmax());
      const auto  hminratio = nodeList.hminratio();
      const auto  nPerh = nodeList.nodesPerSmoothingScale();
      const auto  n = nodeList.numInternalNodes();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        auto& Hi = H(k,i);
        auto& Hideali = Hideal(k,i);
  
        // Complete the zeroth moment
        auto& massZerothMomenti = mZerothMoment(k,i);
        massZerothMomenti = Dimension::rootnu(max(0.0, massZerothMomenti));

        // Determine the current effective number of nodes per smoothing scale.
        const auto currentNodesPerSmoothingScale = (fuzzyEqual(massZerothMomenti, 0.0) ?  // Is this node isolated (no neighbors)?
                                                    0.5*nPerh :
                                                    mWT.equivalentNodesPerSmoothingScale(massZerothMomenti));
        CHECK2(currentNodesPerSmoothingScale > 0.0, "Bad estimate for nPerh effective from kernel: " << currentNodesPerSmoothingScale);

        // The ratio of the desired to current nodes per smoothing scale.
        const auto s = std::min(4.0, std::max(0.25, nPerh/(currentNodesPerSmoothingScale + 1.0e-30)));
        CHECK(s > 0.0);

        // Now determine how to scale the current H to the desired value.
        const auto a = (s < 1.0 ? 
                        0.4*(1.0 + s*s) :
                        0.4*(1.0 + 1.0/(s*s*s)));
        CHECK(1.0 - a + a*s > 0.0);

        // Now a big branch if we're using the normal IdealH or one of the specialized cases.
        if (voronoi) {

          // Find the new normalized target shape
          auto T = mSecondMoment(k,i); // .sqrt();
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
      
          // T *= std::min(4.0, std::max(0.25, 1.0 - a + a*s));
          T *= s;

          // Build the new H tensor
          if (surfacePoint(k, i) == 0) {
            Hideali = (*mHidealFilterPtr)(k, i, Hi, T.Inverse());
          } else {
            Hideali = (*mHidealFilterPtr)(k, i, Hi, Hi);           // Keep the time evolved version for surface points
          }
          Hi = Hideali;     // Since this is the after all our regular state update gotta update the actual H

        } else if (mFixShape) {

          // We're just scaling the fixed H tensor shape, so very close to the normal SPH IdealH algorithm
          Hideali = Hi / (1.0 - a + a*s);

        } else {

          // We scale H in the radial direction only (also force H to be aligned radially).
          CHECK(mRadialOnly);
          const auto nhat = mRadialFunctorPtr->radialUnitVector(k, i, pos(k,i));
          const auto r1 = mRadialFunctorPtr->radialCoordinate(k, i, pos(k,i));
          Hideali = radialEvolution(Hi, nhat, 1.0 - a + a*s, mRadius0(k,i), r1);

        }

        // Apply limiting and set the actual H
        Hideali = (*mHidealFilterPtr)(k, i, Hi, Hideali);
        const auto hev = Hideali.eigenVectors();
        const auto hminEffInv = min(hminInv, max(hmaxInv, hev.eigenValues.minElement())/hminratio);
        Hideali = constructSymTensorWithBoundedDiagonal(hev.eigenValues, hmaxInv, hminEffInv);
        Hideali.rotationalTransform(hev.eigenVectors);
        Hi = Hideali;
      }
    }

  } else {

    // Apply any requested user filtering/alterations to the final H in the case where we're not using the IdealH algorithm
    const auto numNodeLists = dataBase.numFluidNodeLists();
    auto        H = state.fields(HydroFieldNames::H, SymTensor::zero);
    for (auto k = 0u; k < numNodeLists; ++k) {
      const auto& nodeList = H[k]->nodeList();
      const auto  hminInv = safeInvVar(nodeList.hmin());
      const auto  hmaxInv = safeInvVar(nodeList.hmax());
      const auto  hminratio = nodeList.hminratio();
      const auto  n = nodeList.numInternalNodes();
      for (auto i = 0u; i < n; ++i) {
        H(k,i) = (*mHidealFilterPtr)(k, i, H(k,i), H(k,i));
        const auto hev = H(k,i).eigenVectors();
        const auto hminEffInv = min(hminInv, max(hmaxInv, hev.eigenValues.minElement())/hminratio);
        H(k,i) = constructSymTensorWithBoundedDiagonal(hev.eigenValues, hmaxInv, hminEffInv);
        H(k,i).rotationalTransform(hev.eigenVectors);
      }
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

}
