//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- The base class for all ArtificialViscosities in 
// Spheral++.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"

#include "ArtificialViscosity.hh"

#include <algorithm>

using std::vector;
using std::string;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosity<Dimension>::
ArtificialViscosity():
  mClinear(1.0),
  mCquadratic(1.0),
  mQcorrectionOrder(RKOrder::LinearOrder),  
  mBalsaraShearCorrection(false),
  mClMultiplier(FieldStorageType::CopyFields),
  mCqMultiplier(FieldStorageType::CopyFields),
  mShearCorrection(FieldStorageType::CopyFields),
  mCalculateSigma(false),
  mLimiterSwitch(false),
  mEpsilon2(1.0e-2),
  mNegligibleSoundSpeed(1e-10),
  mCsMultiplier(1e-2),
  mEnergyMultiplier(1.0),
  mSigma(FieldStorageType::CopyFields),
  mGradDivVelocity(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosity<Dimension>::
ArtificialViscosity(Scalar Clinear, Scalar Cquadratic, RKOrder QcorrectionOrder):
  mClinear(Clinear),
  mCquadratic(Cquadratic),
  mQcorrectionOrder(QcorrectionOrder), 
  mBalsaraShearCorrection(false),
  mClMultiplier(FieldStorageType::CopyFields),
  mCqMultiplier(FieldStorageType::CopyFields),
  mShearCorrection(FieldStorageType::CopyFields),
  mCalculateSigma(false),
  mLimiterSwitch(false),
  mEpsilon2(1.0e-2),
  mNegligibleSoundSpeed(1e-10),
  mCsMultiplier(1e-2),
  mEnergyMultiplier(1.0),
  mSigma(FieldStorageType::CopyFields),
  mGradDivVelocity(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosity<Dimension>::
~ArtificialViscosity() {
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           const State<Dimension>& state,
           const StateDerivatives<Dimension>& derivs,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
           const typename Dimension::Scalar /*time*/,
           const typename Dimension::Scalar /*dt*/,
           const TableKernel<Dimension>& W) {

  // If needed, calculate grad v and grad div v.
  if (calculateSigma() or limiter()) calculateSigmaAndGradDivV(dataBase,
                                                               state,
                                                               derivs,
                                                               W,
                                                               boundaryBegin,
                                                               boundaryEnd);

  // Start by assuming we're not scaling the linear and quadratic terms by anything.
  // Note we set the last argument here to false, which ensures that if these Fields already exist we don't overwrite
  // their values.  This is needed 'cause some physics options (like the Morris & Monaghan time evolved Q or 
  // Cullen & Dehnen) need to evolve these parameters.
  dataBase.resizeFluidFieldList(mClMultiplier, 1.0, HydroFieldNames::ArtificialViscousClMultiplier, false);
  dataBase.resizeFluidFieldList(mCqMultiplier, 1.0, HydroFieldNames::ArtificialViscousCqMultiplier, false);

  // If we are applying the Balsara shear flow correction term, calculate the
  // per node multiplier.
  dataBase.resizeFluidFieldList(mShearCorrection, 1.0, "Balsara shear correction", true);
  if (this->balsaraShearCorrection()) {

    // State we need.
    const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
    const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);

    // Calculate the shear Q suppression term for all internal fluid nodes.
    const auto& connectivityMap = dataBase.connectivityMap();
    const auto  numNodeLists = connectivityMap.nodeLists().size();
    for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
      const auto ni = mShearCorrection[nodeListi]->numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
        const auto div = fabs(DvDx(nodeListi, i).Trace());
        const auto curl = curlVelocityMagnitude(DvDx(nodeListi, i));
        const auto hmaxinverse = Dimension::rootnu(H(nodeListi, i).Determinant());
        const auto cs = max(negligibleSoundSpeed(), soundSpeed(nodeListi, i));
        CHECK(div >= 0.0);
        CHECK(curl >= 0.0);
        CHECK(hmaxinverse > 0.0);
        CHECK(cs > 0.0);
        mShearCorrection(nodeListi, i) = div/(div + curl + epsilon2()*cs*hmaxinverse);
        CHECK(mShearCorrection(nodeListi, i) >= 0.0 and mShearCorrection(nodeListi, i) <= 1.0);
      }
    }
  }

  // We deliberately do not finalize the boundaries here, depending on the calling environment
  // to know when to do this efficiently.
  for (auto boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(mClMultiplier);
    (*boundItr)->applyFieldListGhostBoundary(mCqMultiplier);
    (*boundItr)->applyFieldListGhostBoundary(mShearCorrection);
  }
}

//------------------------------------------------------------------------------
// Dump the current state of the Q to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  if (calculateSigma()) file.write(mSigma, pathName + "/sigma");
  if (mLimiterSwitch) file.write(gradDivVelocity(), pathName + "/gradDivVelocity");
  file.write(mClMultiplier, pathName + "/ClMultiplier");
  file.write(mCqMultiplier, pathName + "/CqMultiplier");
  file.write(mShearCorrection, pathName + "/shearCorrection");
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  if (calculateSigma()) file.read(mSigma, pathName + "/sigma");
  if (mLimiterSwitch) file.read(mGradDivVelocity, pathName + "/gradDivVelocity");
  file.read(mClMultiplier, pathName + "/ClMultiplier");
  file.read(mCqMultiplier, pathName + "/CqMultiplier");
  file.read(mShearCorrection, pathName + "/shearCorrection");
}

//------------------------------------------------------------------------------
// Figure out the total stress-strain tensor for a given node pair based on 
// the stored average value and the given (position, velocity) pair.
//------------------------------------------------------------------------------
template<typename Dimension>
// inline
typename Dimension::Tensor
ArtificialViscosity<Dimension>::
sigmaij(const typename Dimension::Vector& rij,
        const typename Dimension::Vector& rijUnit,
        const typename Dimension::Vector& vij,
        const typename Dimension::Scalar& hi2,
        const int nodeListID,
        const int nodeID) const {

  REQUIRE(fuzzyEqual(rijUnit.magnitude2(), 1.0));
  REQUIRE(distinctlyGreaterThan(hi2, 0.0));
  REQUIRE(nodeListID >= 0 and nodeListID < (int)mSigma.size());
  REQUIRE(nodeID >= 0 and nodeID < (int)mSigma[nodeListID]->nodeListPtr()->numNodes());

  // Get the rotational transformations.
  const auto R = rotationMatrix(rijUnit);  // Gets -1
  const auto Rinverse = R.Transpose();

  // Rotate the velocity vector into the rotated frame.  We'll label this 
  // as deltaSigma, since that's what we're going to construct with this
  // rotated velocity.
  auto deltaSigma = R*vij;

  // Construct the pairwise column of the sigma tensor as a vector.
  const auto xprime2 = rij.magnitude2();
  const auto norm = sqrt(xprime2)/(xprime2 + mEpsilon2*hi2);
  deltaSigma *= norm;

  // Now combine the background and pairwise values in the rotated frame,
  // then rotate back to the lab frame and return the result.
  auto sigma = (*mSigma[nodeListID])(nodeID);
  sigma.rotationalTransform(R);
  sigma.setColumn(0, deltaSigma);
  sigma.rotationalTransform(Rinverse);
  return sigma;
}

//------------------------------------------------------------------------------
// Compute the internal background sigma and grad-div-v fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ArtificialViscosity<Dimension>::
calculateSigmaAndGradDivV(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& /*derivs*/,
                          const TableKernel<Dimension>& W,
                          typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
                          typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd) {

  const double tiny = 1.0e-10;

  // Verify that the internal FieldLists mSigma and mGradDivVelocity are properly
  // sized.
  dataBase.resizeFluidFieldList(mSigma, Tensor::zero, "sigmaQ");
  dataBase.resizeFluidFieldList(mGradDivVelocity, Vector::zero, "gradDivVelocity");

  // Zero out the gradient and grad div v FieldLists.
  mSigma.Zero();
  mGradDivVelocity.Zero();

  // Get the necessary state.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto Hfield = state.fields(HydroFieldNames::H, SymTensor::zero);

  // Prepare FieldLists to accumulate the normalizations in.
  FieldList<Dimension, Tensor> sigNormalization = dataBase.newFluidFieldList(Tensor::zero, "sigma normalization");
  FieldList<Dimension, Scalar> gdvNormalization = dataBase.newFluidFieldList(0.0, "grad div v normalization");

  // Grab the connectivity map from the DataBase.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = dataBase.numFluidNodeLists();
  CONTRACT_VAR(nodeLists);
  CHECK(nodeLists.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, Wj, gWj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto sigma_thread = mSigma.threadCopy(threadStack);
    auto gradDivVelocity_thread = mGradDivVelocity.threadCopy(threadStack);
    auto sigNormalization_thread = sigNormalization.threadCopy(threadStack);
    auto gdvNormalization_thread = gdvNormalization.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // State for node i.
      const auto& ri = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& Hi = Hfield(nodeListi, i);
      const auto  weighti = mass(nodeListi, i)/rho(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      // The derivatives for i.
      auto& sigmai = sigma_thread(nodeListi, i);
      auto& gradDivVelocityi = gradDivVelocity_thread(nodeListi, i);
      auto& sigNormalizationi = sigNormalization_thread(nodeListi, i);
      auto& gdvNormalizationi = gdvNormalization_thread(nodeListi, i);

      // State for node j.
      const auto& rj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& Hj = Hfield(nodeListj, j);
      const auto  weightj = mass(nodeListj, j)/rho(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();

      // The derivatives for j.
      auto& sigmaj = sigma_thread(nodeListj, j);
      auto& gradDivVelocityj = gradDivVelocity_thread(nodeListj, j);
      auto& sigNormalizationj = sigNormalization_thread(nodeListj, j);
      auto& gdvNormalizationj = gdvNormalization_thread(nodeListj, j);

      // Get an estimate of the smoothing scale for i and j.
      CHECK(distinctlyGreaterThan(Hdeti, 0.0));
      CHECK(distinctlyGreaterThan(Hdetj, 0.0));
      const auto frach = 0.02/(Dimension::rootnu(Hdeti) + Dimension::rootnu(Hdetj));
      const auto frach2 = frach*frach;
      CHECK(distinctlyGreaterThan(frach, 0.0));
      CHECK(distinctlyGreaterThan(frach2, 0.0));

      // Compute the relative distance and eta for i & j.
      const auto rij = ri - rj;
      const auto etai = Hi*rij;
      const auto etamagi = etai.magnitude();
      const auto Hetai = Hi*etai.unitVector();
      std::tie(Wi, gWi) = W.kernelAndGradValue(etamagi, Hdeti);

      const auto etaj = Hj*rij;
      const auto etamagj = etaj.magnitude();
      const auto Hetaj = Hj*etaj.unitVector();
      std::tie(Wj, gWj) = W.kernelAndGradValue(etamagj, Hdetj);

      // Sum this pairs contribution to grad v for both i & j.
      const auto vij = vi - vj;

      const auto wij = weightj*Wi;
      const auto gwij = weightj*Hetai*gWi;

      const auto wji = weighti*Wj;
      const auto gwji = -weighti*Hetaj*gWj;

      const auto rjiUnit = -rij.unitVector();
      const auto R = rotationMatrix(rjiUnit);
      const auto Rinverse = R.Transpose();
      const auto rij2 = rij.magnitude2();
      const auto rijmag = sqrt(rij2);

      const auto dxp = rijmag/(rij2 + frach2);
      const auto dfdxpcol = -dxp*(R*vij);
      Tensor dfdxp;
      dfdxp.setColumn(0, dfdxpcol);
      dfdxp.rotationalTransform(Rinverse);

      sigmai += wij*dfdxp;
      sigmaj += wji*dfdxp;
            
      // Update the sigma normalization.
      auto tweighti = constructTensorWithColumnValue<Tensor, 0>(wij);
      auto tweightj = constructTensorWithColumnValue<Tensor, 0>(wji);
      tweighti.rotationalTransform(Rinverse);
      tweightj.rotationalTransform(-Rinverse);
      inPlaceAbsAdd<Tensor>(sigNormalizationi, tweighti);
      inPlaceAbsAdd<Tensor>(sigNormalizationj, tweightj);

      // Sum this pairs contribution to the second derivative of v.
      const auto thpt = vij.dot(rij)/(rij2 + frach2);

      gradDivVelocityi += thpt*gwij;
      gdvNormalizationi += wij;

      gradDivVelocityj += thpt*gwji;
      gdvNormalizationj += wji;
    }

    // Reduce the thread values to the masters.
    threadReduceFieldLists<Dimension>(threadStack);

  } // OpenMP parallel region

  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = mSigma[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // The derivatives for i.
      auto& sigmai = mSigma(nodeListi, i);
      auto& gradDivVelocityi = mGradDivVelocity(nodeListi, i);
      auto& sigNormalizationi = sigNormalization(nodeListi, i);
      auto& gdvNormalizationi = gdvNormalization(nodeListi, i);

      // Complete the sigma calculation for i.
      CHECK(*std::min_element(sigNormalizationi.begin(), sigNormalizationi.end()) >= 0.0);
      tensorElementWiseDivide<Tensor>(sigmai, sigNormalizationi + tiny*Tensor::one);

      // Now limit to just negative eigen-values.  This is 'cause we only
      // care about convergent geometries for the Q.
      const auto sigmai_s = sigmai.Symmetric();
      const auto sigmai_a = sigmai.SkewSymmetric();
      auto eigeni = sigmai_s.eigenVectors();
      sigmai = constructTensorWithMinDiagonal(eigeni.eigenValues, 0.0);
      sigmai.rotationalTransform(eigeni.eigenVectors);
      sigmai += sigmai_a;

      // Complete grad div v.
      CHECK(gdvNormalizationi >= 0.0);
      gradDivVelocityi /= gdvNormalizationi + tiny;
    }
  }

  // Apply boundary conditions
  for (auto boundItr = boundaryBegin;
       boundItr < boundaryEnd;
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(mSigma);
    (*boundItr)->applyFieldListGhostBoundary(mGradDivVelocity);
  }
  // for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
  //      boundItr != boundaryEnd;
  //      ++boundItr) {
  //   (*boundItr)->finalizeGhostBoundary();
  // }

}

}
