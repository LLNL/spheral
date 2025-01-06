//---------------------------------Spheral++----------------------------------//
// ASPHClassicSmoothingScale
//
// Implements our classic ASPH algorithm (2010 SPHERIC proceedings style)
//
// Created by JMO, Fri Jan  3 10:48:43 PST 2025
//----------------------------------------------------------------------------//
#include "SmoothingScale/ASPHClassicSmoothingScale.hh"
#include "SmoothingScale/SmoothingScaleUtilities.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/Timer.hh"

#include <cmath>
#include <vector>

namespace Spheral {

using std::min;
using std::max;
using std::abs;
using std::vector;

namespace {

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ASPHClassicSmoothingScale<Dimension>::
ASPHClassicSmoothingScale(const HEvolutionType HUpdate,
                          const TableKernel<Dimension>& W):
  SmoothingScaleBase<Dimension>(HUpdate),
  mWT(W),
  mZerothMoment(FieldStorageType::CopyFields),
  mFirstMoment(FieldStorageType::CopyFields),
  mSecondMoment(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHClassicSmoothingScale<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  // Make sure our FieldLists are correctly sized.
  SmoothingScaleBase<Dimension>::initializeProblemStartup(dataBase);
  dataBase.resizeFluidFieldList(mZerothMoment, 0.0, HydroFieldNames::massZerothMoment, false);
  dataBase.resizeFluidFieldList(mFirstMoment, Vector::zero, HydroFieldNames::massFirstMoment, false);
  dataBase.resizeFluidFieldList(mSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
}

//------------------------------------------------------------------------------
// Register derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHClassicSmoothingScale<Dimension>::
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
ASPHClassicSmoothingScale<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("ASPHClassicSmoothingScaleDerivs");

  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  etaMax = mWT.kernelExtent();

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
  auto  massSecondMoment = derivs.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(massZerothMoment.size() == numNodeLists);
  CHECK(massFirstMoment.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar mi, mj, ri, rj, mRZi, mRZj, rhoi, rhoj, WSPHi, WSPHj, etaMagi, etaMagj, fweightij, fispherical, fjspherical;
    Vector xij, etai, etaj;
    SymTensor xijdyad;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto massZerothMoment_thread = massZerothMoment.threadCopy(threadStack);
    auto massFirstMoment_thread = massFirstMoment.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      mi = mass(nodeListi, i);
      rhoi = massDensity(nodeListi, i);
      const auto& xi = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);

      auto& massZerothMomenti = massZerothMoment_thread(nodeListi, i);
      auto& massFirstMomenti = massFirstMoment_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);

      // Get the state for node j
      mj = mass(nodeListj, j);
      rhoj = massDensity(nodeListj, j);
      const auto& xj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);

      auto& massZerothMomentj = massZerothMoment_thread(nodeListj, j);
      auto& massFirstMomentj = massFirstMoment_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);

      // Node displacement.
      xij = xi - xj;
      etai = Hi*xij;
      etaj = Hj*xij;
      etaMagi = etai.magnitude();
      etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Compute the node-node weighting
      fweightij = 1.0;
      fispherical = 1.0;
      fjspherical = 1.0;
      if (nodeListi != nodeListj) {
        if (GeometryRegistrar::coords() == CoordinateType::RZ) {
          ri = abs(xi.y());
          rj = abs(xj.y());
          mRZi = mi/(2.0*M_PI*ri);
          mRZj = mj/(2.0*M_PI*rj);
          fweightij = mRZj*rhoi/(mRZi*rhoj);
        } else {
          fweightij = mj*rhoi/(mi*rhoj);
        }
      } else if (GeometryRegistrar::coords() == CoordinateType::Spherical) {
        const auto eii = Hi.xx()*xi.x();
        const auto eji = Hi.xx()*xj.x();
        const auto ejj = Hj.xx()*xj.x();
        const auto eij = Hj.xx()*xi.x();
        fispherical = (eii > etaMax ? 1.0 :
                       eii < eji ? 2.0 :
                       0.0);
        fjspherical = (ejj > etaMax ? 1.0 :
                       ejj < eij ? 2.0 :
                       0.0);
      }

      // Symmetrized kernel weight
      WSPHi = mWT.kernelValueSPH(etaMagi);
      WSPHj = mWT.kernelValueSPH(etaMagj);

      // Moments of the node distribution -- used for the ideal H calculation.
      massZerothMomenti +=     fweightij*WSPHi * fispherical;
      massZerothMomentj += 1.0/fweightij*WSPHj * fjspherical;
      massFirstMomenti -=     fweightij*WSPHi*etai;
      massFirstMomentj += 1.0/fweightij*WSPHj*etaj;
      // massSecondMomenti +=     fweightij*WSPHi*WSPHi*etai.unitVector().selfdyad();
      // massSecondMomentj += 1.0/fweightij*WSPHj*WSPHj*etaj.unitVector().selfdyad();
      xijdyad = xij.selfdyad()*safeInvVar(FastMath::pow5(xij.magnitude()));
      massSecondMomenti +=     fweightij*WSPHi*WSPHi*xijdyad;
      massSecondMomentj += 1.0/fweightij*WSPHj*WSPHj*xijdyad;
    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region

  // Finish up the derivatives now that we've walked all pairs
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hminInv = safeInvVar(nodeList.hmin());
    const auto  hmaxInv = safeInvVar(nodeList.hmax());
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& Hi = H(nodeListi, i);
      const auto& DvDxi = DvDx(nodeListi, i);
      auto&       massZerothMomenti = massZerothMoment(nodeListi, i);
      const auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto&       DHDti = DHDt(nodeListi, i);
      auto&       Hideali = Hideal(nodeListi, i);

      // Time derivative of H
      DHDti = SmoothingScaleDetail::smoothingScaleDerivative(Hi, DvDxi);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      massZerothMomenti = Dimension::rootnu(max(0.0, massZerothMomenti));

      // Determine the current effective number of nodes per smoothing scale.
      const auto currentNodesPerSmoothingScale = (fuzzyEqual(massZerothMomenti, 0.0) ?  // Is this node isolated (no neighbors)?
                                                  0.5*nPerh :
                                                  mWT.equivalentNodesPerSmoothingScale(massZerothMomenti));
      CHECK2(currentNodesPerSmoothingScale > 0.0, "Bad estimate for nPerh effective from kernel: " << currentNodesPerSmoothingScale);

      // The (limited) ratio of the current to desired nodes per smoothing scale.
      const Scalar s = min(4.0, max(0.25, currentNodesPerSmoothingScale/nPerh));
      CHECK(s > 0.0);

      // Start with the sqrt of the second moment in eta space
      Hideali = massSecondMomenti.sqrt();
      auto eigenT = Hideali.eigenVectors();

      // Ensure we don't have any degeneracies (zero eigen values)
      const auto Tmax = std::max(1.0, eigenT.eigenValues.maxElement());
      auto fscale = 1.0;
      for (auto k = 0u; k < Dimension::nDim; ++k) {
        eigenT.eigenValues[k] = max(eigenT.eigenValues[k], 0.01*Tmax);
        fscale *= eigenT.eigenValues[k];
      }
      CHECK(fscale > 0.0);

      // Compute the scaling to get us closer to the target n per h, and build the transformation tensor
      fscale = 1.0/Dimension::rootnu(fscale);    // inverse length, same as H!
      eigenT.eigenValues *= fscale;
      Hideali = constructSymTensorWithDiagonal(eigenT.eigenValues);
      Hideali.rotationalTransform(eigenT.eigenVectors);
      CHECK(fuzzyEqual(Hideali.Determinant(), 1.0, 1.0e-8));

      // Initial vote for Hideal
      Hideali *= s*Dimension::rootnu(Hi.Determinant());
      // Hideali = (T*Hi).Symmetric();

      // Apply limiting
      const auto hev = Hideali.eigenVectors();
      const auto hminEffInv = min(hminInv, max(hmaxInv, hev.eigenValues.minElement())/hminratio);
      Hideali = constructSymTensorWithBoundedDiagonal(hev.eigenValues, hmaxInv, hminEffInv);
      Hideali.rotationalTransform(hev.eigenVectors);
    }
  }
  TIME_END("ASPHClassicSmoothingScaleDerivs");
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ASPHClassicSmoothingScale<Dimension>::
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
ASPHClassicSmoothingScale<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  SmoothingScaleBase<Dimension>::restoreState(file, pathName);
  file.read(mZerothMoment, pathName + "/zerothMoment");
  file.read(mFirstMoment, pathName + "/firstMoment");
  file.read(mSecondMoment, pathName + "/secondMoment");
}

}
