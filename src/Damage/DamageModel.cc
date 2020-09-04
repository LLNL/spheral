//---------------------------------Spheral++----------------------------------//
// DamageModel -- Base class for the damage physics models.
// This class just provides the basic interface for damage models, and does 
// not fill out the complete physics package interface.
//
// Created by JMO, Thu Sep 29 13:31:57 PDT 2005
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DamageModel.hh"
#include "YoungsModulusPolicy.hh"
#include "LongitudinalSoundSpeedPolicy.hh"
#include "DamagedSoundSpeedPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/safeInv.hh"

#include "boost/shared_ptr.hpp"

#include <string>
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
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DamageModel<Dimension>::
DamageModel(SolidNodeList<Dimension>& nodeList,
            const TableKernel<Dimension>& W,
            const double crackGrowthMultiplier,
            const EffectiveFlawAlgorithm flawAlgorithm,
            const FlawStorageType& flaws):
  Physics<Dimension>(),
  mFlaws(SolidFieldNames::flaws, flaws),
  mEffectiveFlaws(SolidFieldNames::effectiveFlaws, nodeList),
  mNodeList(nodeList),
  mW(W),
  mCrackGrowthMultiplier(crackGrowthMultiplier),
  mEffectiveFlawAlgorithm(flawAlgorithm),
  mYoungsModulus(SolidFieldNames::YoungsModulus, nodeList),
  mLongitudinalSoundSpeed(SolidFieldNames::longitudinalSoundSpeed, nodeList),
  mExcludeNode("Nodes excluded from damage", nodeList, 0),
  mCriticalNodesPerSmoothingScale(0.99),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DamageModel<Dimension>::
~DamageModel() {
}

//------------------------------------------------------------------------------
// Compute the generic Grady-Kipp (ala Benz-Asphaug) scalar damage time 
// derivative.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
computeScalarDDDt(const DataBase<Dimension>& /*dataBase*/,
                  const State<Dimension>& state,
                  const Scalar /*time*/,
                  const Scalar /*dt*/,
                  Field<Dimension, Scalar>& DDDt) const {

  // Pre-conditions.
  REQUIRE(DDDt.nodeListPtr() == &mNodeList);
  REQUIRE(mFlaws.nodeListPtr() == &mNodeList);

  // Get the state fields.
  const auto  clKey = State<Dimension>::buildFieldKey(SolidFieldNames::longitudinalSoundSpeed, mNodeList.name());
  const auto  HKey = State<Dimension>::buildFieldKey(HydroFieldNames::H, mNodeList.name());
  const auto& cl = state.field(clKey, 0.0);
  const auto& H = state.field(HKey, SymTensor::zero);

  // Constant multiplicative parameter for the crack growth.
  const auto A = mCrackGrowthMultiplier / mW.kernelExtent();

  // Iterate over the internal nodes.
  const auto ni = mNodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < ni; ++i) {
    if (mExcludeNode(i) == 1) {

      DDDt(i) = 0.0;

    } else {

      const double hrInverse = Dimension::rootnu(H(i).Determinant());
      DDDt(i) = A * cl(i) * hrInverse;

    }
    CHECK(DDDt(i) >= 0.0);
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  for (int i = 0; i != mNodeList.numInternalNodes(); ++i) {
    if (mExcludeNode(i) == 1) {
      ENSURE(DDDt(i) == 0.0);
    }
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // Register Youngs modulus, the longitudinal sound speed, and the effective flaw strains.
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;
  PolicyPointer EPolicy(new YoungsModulusPolicy<Dimension>());
  PolicyPointer clPolicy(new LongitudinalSoundSpeedPolicy<Dimension>());
  state.enroll(mYoungsModulus, EPolicy);
  state.enroll(mLongitudinalSoundSpeed, clPolicy);
  state.enroll(mEffectiveFlaws);

  // Set the initial values for the Youngs modulus, sound speed, pressure, and effective flaws.
  typename StateDerivatives<Dimension>::PackageList dummyPackages;
  StateDerivatives<Dimension> derivs(dataBase, dummyPackages);
  EPolicy->update(state.key(mYoungsModulus), state, derivs, 1.0, 0.0, 0.0);
  clPolicy->update(state.key(mLongitudinalSoundSpeed), state, derivs, 1.0, 0.0, 0.0);
}

//------------------------------------------------------------------------------
// preStepInitialize -- for this class this means determine the effective flaw
// failure strains.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {

  // If we're just using the "FullSpectrum" of flaws, we don't reduce the 
  // flaws for each node to a single value at all.  In that case there's no
  // work to do.
  if (mEffectiveFlawAlgorithm != EffectiveFlawAlgorithm::FullSpectrumFlaws) {

    // Get the state fields.
    const auto  posKey = State<Dimension>::buildFieldKey(HydroFieldNames::position, mNodeList.name());
    const auto  HKey = State<Dimension>::buildFieldKey(HydroFieldNames::H, mNodeList.name());
    const auto  flawKey = State<Dimension>::buildFieldKey(SolidFieldNames::effectiveFlaws, mNodeList.name());
    const auto& positions = state.field(posKey, Vector::zero);
    const auto& H = state.field(HKey, SymTensor::zero);
    auto&       effectiveFlaws = state.field(flawKey, 0.0);

    // Iterate over the nodes and compute our effective flaws!
    const auto ni = mNodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto& flawsi = mFlaws(i);
      if (flawsi.size() > 0) {

        switch(mEffectiveFlawAlgorithm) {

        case EffectiveFlawAlgorithm::MinFlaw:
          effectiveFlaws(i) = *min_element(flawsi.begin(), flawsi.end());
          break;

        case EffectiveFlawAlgorithm::MaxFlaw:
          effectiveFlaws(i) = *max_element(flawsi.begin(), flawsi.end());
          break;

        case EffectiveFlawAlgorithm::InverseSumFlaws:
        case EffectiveFlawAlgorithm::SampledFlaws:
          effectiveFlaws(i) = 0.0;
          for (auto itr = flawsi.begin(); itr != flawsi.end(); ++itr) effectiveFlaws(i) += 1.0/(*itr);
          effectiveFlaws(i) = flawsi.size()/effectiveFlaws(i);
          break;

        default:
          CHECK(false);
          break;

        }

      }
    }

    // For the sampled flaws case, we have to communicate the local sums and then
    // do the weighted interpolation.
    if (mEffectiveFlawAlgorithm == EffectiveFlawAlgorithm::SampledFlaws) {
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto& nodeLists = connectivityMap.nodeLists();
      const auto  numNodeLists = nodeLists.size();
      const auto  nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), &mNodeList));
      CONTRACT_VAR(nodeListi);
      CONTRACT_VAR(numNodeLists);
      CHECK(nodeListi < numNodeLists);

      // Invert the flaws again, and remove the number of flaws normalization.
#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
        effectiveFlaws(i) = mFlaws(i).size()/effectiveFlaws(i);
      }

      // Apply ghost boundaries.
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyGhostBoundary(effectiveFlaws);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

      // Make a copy of the inverse sum for each node.
      Field<Dimension, Scalar> flawCopy("flaws copy", effectiveFlaws);

      // Go over everybody again.
#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {

        // The self-contribution.
        const auto& ri = positions(i);
        const auto& Hi = H(i);
        const auto  Hdeti = Hi.Determinant();
        const auto  Wi = mW(0.0, Hdeti);
        effectiveFlaws(i) *= Wi;
        auto normalization = Wi;

        // Now everyone else's contribution.
        const auto& connectivity = connectivityMap.connectivityForNode(&mNodeList, i)[nodeListi];
        for (auto jItr = connectivity.begin(); jItr < connectivity.end(); ++jItr) {
          const int j = *jItr;
          CHECK(j < mNodeList.numNodes());
          const auto rij = ri - positions(j);
          const auto Wi = mW(Hi*rij, Hdeti);
          effectiveFlaws(i) += Wi*flawCopy(j);
          normalization += Wi;
        }

        // Normalize and we're done with this node.
        effectiveFlaws(i) = normalization/effectiveFlaws(i);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Post-state update jobs.
//------------------------------------------------------------------------------
template<typename Dimension>
void 
DamageModel<Dimension>::
postStateUpdate(const Scalar /*time*/, 
                const Scalar /*dt*/,
                const DataBase<Dimension>& /*dataBase*/, 
                State<Dimension>& /*state*/,
                StateDerivatives<Dimension>& /*derivatives*/) {

//   typedef typename SymTensor::EigenStructType EigenStruct;

//   // Connectivity info.
//   const auto& connectivityMap = dataBase.connectivityMap();
//   const auto& nodeLists = connectivityMap.nodeLists();
//   const auto  numNodeLists = nodeLists.size();
//   const auto  nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), &mNodeList));
//   CHECK(nodeListi < numNodeLists);

//   // Grab the state.
//   typedef typename State<Dimension>::KeyType Key;
//   const auto& x = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::position, mNodeList.name()), Vector::zero);
//   const auto& D = state.field(State<Dimension>::buildFieldKey(SolidFieldNames::tensorDamage, mNodeList.name()), SymTensor::zero);
//   const auto& m = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::mass, mNodeList.name()), 0.0);
//   const auto& rho = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, mNodeList.name()), 0.0);
//   auto&       H = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::H, mNodeList.name()), SymTensor::zero);
//   auto&       S = state.field(State<Dimension>::buildFieldKey(SolidFieldNames::deviatoricStress, mNodeList.name()), SymTensor::zero);

//   // Force the deviatoric stress of the damaged points to limit to the local average.
//   Field<Dimension, SymTensor> Snew("S sampled", mNodeList);
// #pragma omp parallel for
//   for (auto i = 0; i < mNodeList.numInternalNodes(); ++i) {
//     const auto& xi = x(i);
//     const auto& Hi = H(i);
//     const auto  Hdeti = Hi.Determinant();
//     const auto  Dmax = max(0.0, min(1.0, D(i).eigenValues().maxElement()));
//     const auto& connectivity = connectivityMap.connectivityForNode(&mNodeList, i)[nodeListi];
//     if (connectivity.size() > 0 and Dmax > 0.0) {
//       const auto weighti = m(i)*safeInv(rho(i));
//       auto       weightSum = (1.0 - Dmax)*weighti*mW(0.0, Hdeti);
//       Snew(i) = weightSum*S(i);
//       for (auto jItr = connectivity.begin(); jItr < connectivity.end(); ++jItr) {
//         const auto j = *jItr;
//         const auto weightj = m(j)*safeInv(rho(j));
//         const auto xij = xi - x(j);
//         const auto wj = max(0.0, min(1.0, 1.0 - D(j).eigenValues().maxElement())) * weightj * mW(Hi*xij, Hi);
//         weightSum += wj;
//         Snew(i) += wj*S(j);
//       }
//       Snew(i) *= safeInv(weightSum, 1.0e-10);
//       Snew(i) = (1.0 - Dmax)*S(i) + Dmax*Snew(i);
//     } else {
//       Snew(i) = S(i);
//     }
//   }

//   // We require the H tensors of nodes approach unit aspect ratio as the damage
//   // is increased.
//   // Iterate over our nodes and limit H's as need be.
// #pragma omp parallel for
//   for (auto i = 0; i < mNodeList.numInternalNodes(); ++i) {
//     const auto& Di = D(i);
//     auto&       Hi = H(i);
//     const auto  Dmax = max(0.0, min(1.0, Di.eigenValues().maxElement()));
//     if (Dmax > mNodeList.hminratio()) {
//       const auto Heigen = Hi.Inverse().eigenVectors();
//       const auto hmineff = Dmax*Heigen.eigenValues.maxElement();
//       auto       HnewInv = constructSymTensorWithMaxDiagonal(Heigen.eigenValues, hmineff);
//       HnewInv.rotationalTransform(Heigen.eigenVectors);
//       Hi = HnewInv.Inverse();
//     }
//   }

//   // Apply boundary conditions since we've changed the state.
//   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
//        boundaryItr != this->boundaryEnd();
//        ++boundaryItr) {
//     (*boundaryItr)->applyGhostBoundary(S);
//     (*boundaryItr)->applyGhostBoundary(H);
//   }

}

//------------------------------------------------------------------------------
// Cull the flaws on each node to the single weakest one.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
cullToWeakestFlaws() {
#pragma omp parallel for
  for (auto i = 0u; i < mNodeList.numInternalNodes(); ++i) {
    auto& flaws = mFlaws[i];
    if (flaws.size() > 0) {
      const auto maxVal = *max_element(flaws.begin(), flaws.end());
      flaws = vector<double>(maxVal);
    }
  }
}

//------------------------------------------------------------------------------
// Access the set of nodes to be excluded from damage.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<int>
DamageModel<Dimension>::
excludeNodes() const {
  vector<int> result;
  for (auto i = 0u; i != mNodeList.numInternalNodes(); ++i) {
    if (mExcludeNode(i) == 1.0) result.push_back(i);
  }
  return result;
}

template<typename Dimension>
void
DamageModel<Dimension>::
excludeNodes(vector<int> ids) {
  mExcludeNode = 0;
  for (vector<int>::const_iterator itr = ids.begin();
       itr != ids.end();
       ++itr) {
    REQUIRE(*itr >= 0 && *itr < mNodeList.numInternalNodes());
    mExcludeNode(*itr) = 1;
  }
}

//------------------------------------------------------------------------------
// Compute a Field with the sum of the activation energies per node.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
DamageModel<Dimension>::
sumActivationEnergiesPerNode() const {
  Field<Dimension, Scalar> result("Sum activation energies", mNodeList);
  for (auto i = 0u; i != mNodeList.numInternalNodes(); ++i) {
    const vector<double>& flaws = mFlaws(i);
    for (vector<double>::const_iterator flawItr = flaws.begin();
         flawItr != flaws.end();
         ++flawItr) {
      result(i) += *flawItr;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Compute a Field with the number of flaws per node.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
DamageModel<Dimension>::
numFlawsPerNode() const {
  Field<Dimension, Scalar> result("num flaws", mNodeList);
  for (auto i = 0u; i != mNodeList.numInternalNodes(); ++i) {
    result(i) = flawsForNode(i).size();
  }
  return result;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mCrackGrowthMultiplier, pathName + "/crackGrowthMultiplier");
  file.write(mFlaws, pathName + "/flaws");
  file.write(mExcludeNode, pathName + "/excludeNode");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamageModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mCrackGrowthMultiplier, pathName + "/crackGrowthMultiplier");
  file.read(mFlaws, pathName + "/flaws");
  file.read(mExcludeNode, pathName + "/excludeNode");
}

}

