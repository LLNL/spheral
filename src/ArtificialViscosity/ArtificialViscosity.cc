//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- The base class for all ArtificialViscosities in 
// Spheral++.
//----------------------------------------------------------------------------//

#include <algorithm>

#include "ArtificialViscosity.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"
#include "FileIO/FileIO.hh"

using namespace std;

namespace Spheral {
namespace ArtificialViscositySpace {

using std::vector;
using NodeSpace::NodeList;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using BoundarySpace::Boundary;
using KernelSpace::TableKernel;
using NeighborSpace::Neighbor;
using FileIOSpace::FileIO;
using NeighborSpace::ConnectivityMap;
using NodeSpace::FluidNodeList;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosity<Dimension>::
ArtificialViscosity():
  mClinear(1.0),
  mCquadratic(1.0),
  mBalsaraShearCorrection(false),
  mShearMultiplier(FieldSpace::Copy),
  mReducingViscosityCorrection(false),
  mReducingViscosityMultiplierQ(FieldSpace::Copy),
  mReducingViscosityMultiplierL(FieldSpace::Copy),
  mCalculateSigma(false),
  mLimiterSwitch(false),
  mEpsilon2(1.0e-2),
  mNegligibleSoundSpeed(1e-10),
  mCsMultiplier(1e-2),
  mEnergyMultiplier(1.0),
  mSigma(FieldSpace::Copy),
  mGradDivVelocity(FieldSpace::Copy),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
ArtificialViscosity<Dimension>::
ArtificialViscosity(Scalar Clinear, Scalar Cquadratic):
  mClinear(Clinear),
  mCquadratic(Cquadratic),
  mBalsaraShearCorrection(false),
  mShearMultiplier(FieldSpace::Copy),
  mReducingViscosityCorrection(false),
  mReducingViscosityMultiplierQ(FieldSpace::Copy),
  mReducingViscosityMultiplierL(FieldSpace::Copy),
  mCalculateSigma(false),
  mLimiterSwitch(false),
  mEpsilon2(1.0e-4),
  mNegligibleSoundSpeed(1e-10),
  mCsMultiplier(1e-2),
  mEnergyMultiplier(1.0),
  mSigma(FieldSpace::Copy),
  mGradDivVelocity(FieldSpace::Copy),
  mRestart(DataOutput::registerWithRestart(*this)) {
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
           const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

    // If needed, calculate grad v and grad div v.
    if (calculateSigma() or limiter()) calculateSigmaAndGradDivV(dataBase,
                                                               state,
                                                               derivs,
                                                               W,
                                                               boundaryBegin,
                                                               boundaryEnd);

    // add rvAlpha to the FluidFieldList if it doesn't already exist
    if (reducingViscosityCorrection()) {
      dataBase.resizeFluidFieldList(mReducingViscosityMultiplierQ, 0.1, HydroFieldNames::reducingViscosityMultiplierQ, false);
      dataBase.resizeFluidFieldList(mReducingViscosityMultiplierL, 0.1, HydroFieldNames::reducingViscosityMultiplierL, false);
    } else {
      dataBase.resizeFluidFieldList(mReducingViscosityMultiplierQ, 1.0, HydroFieldNames::reducingViscosityMultiplierQ, false);
      dataBase.resizeFluidFieldList(mReducingViscosityMultiplierL, 1.0, HydroFieldNames::reducingViscosityMultiplierL, false);
    }
    for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
         boundItr != boundaryEnd;
         ++boundItr) {
      (*boundItr)->applyFieldListGhostBoundary(mReducingViscosityMultiplierQ);
      (*boundItr)->applyFieldListGhostBoundary(mReducingViscosityMultiplierL);
    }
    
    // If we are applying the Balsara shear flow correction term, calculate the
    // per node multiplier.
    if (balsaraShearCorrection()) {

      // Begin by verifying that the internal FieldList mQultiplier is correctly
      // sized.
      dataBase.resizeFluidFieldList(mShearMultiplier, 0.0, "sigmaShear");

      // State we need.
      const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
      const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
      const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);

      // Calculate the shear Q suppression term for all internal fluid nodes.
      const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
      const int numNodeLists = connectivityMap.nodeLists().size();
      for (int nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
        for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
             iItr != connectivityMap.end(nodeListi);
             ++iItr) {
          const int i = *iItr;
          const Scalar div = fabs(DvDx(nodeListi, i).Trace());
          const Scalar curl = curlVelocityMagnitude(DvDx(nodeListi, i));
          const Scalar hmaxinverse = Dimension::rootnu(H(nodeListi, i).Determinant());
          const Scalar cs = max(negligibleSoundSpeed(), soundSpeed(nodeListi, i));
          CHECK(div >= 0.0);
          CHECK(curl >= 0.0);
          CHECK(hmaxinverse > 0.0);
          CHECK(cs > 0.0);
          mShearMultiplier(nodeListi, i) = div/(div + curl + epsilon2()*cs*hmaxinverse);
          CHECK(mShearMultiplier(nodeListi, i) >= 0.0 and mShearMultiplier(nodeListi, i) <= 1.0);
        }
      }

      // Apply boundary conditions to fill out the shear corrections on ghost nodes.
      // We deliberately do not finalize the boundaries here, depending on the calling environment
      // to know when to do this efficiently.
      for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
           boundItr != boundaryEnd;
           ++boundItr) {
        (*boundItr)->applyFieldListGhostBoundary(mShearMultiplier);
      }
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
  if (mBalsaraShearCorrection) file.write(mShearMultiplier, pathName + "/shearMultiplier");
  file.write(mReducingViscosityMultiplierQ, pathName + "/reducingViscosityMultiplierQ");
  file.write(mReducingViscosityMultiplierL, pathName + "/reducingViscosityMultiplierL");
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
  if (mBalsaraShearCorrection) file.read(mShearMultiplier, pathName + "/shearMultiplier");
  file.read(mReducingViscosityMultiplierQ, pathName + "/reducingViscosityMultiplierQ");
  file.read(mReducingViscosityMultiplierL, pathName + "/reducingViscosityMultiplierL");
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
  REQUIRE(nodeListID >= 0 and nodeListID < mSigma.size());
  REQUIRE(nodeID >= 0 and nodeID < mSigma[nodeListID]->nodeListPtr()->numNodes());

  // Get the rotational transformations.
  const Tensor R = rotationMatrix(rijUnit);  // Gets -1
  const Tensor Rinverse = R.Transpose();

  // Rotate the velocity vector into the rotated frame.  We'll label this 
  // as deltaSigma, since that's what we're going to construct with this
  // rotated velocity.
  Vector deltaSigma = R*vij;

  // Construct the pairwise column of the sigma tensor as a vector.
  const Scalar xprime2 = rij.magnitude2();
  const Scalar norm = sqrt(xprime2)/(xprime2 + mEpsilon2*hi2);
  deltaSigma *= norm;

  // Now combine the background and pairwise values in the rotated frame,
  // then rotate back to the lab frame and return the result.
  Tensor sigma = (*mSigma[nodeListID])(nodeID);
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
                          const StateDerivatives<Dimension>& derivs,
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
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, SymTensor> Hfield = state.fields(HydroFieldNames::H, SymTensor::zero);

  // Prepare FieldLists to accumulate the normalizations in.
  FieldList<Dimension, Tensor> sigNormalization = dataBase.newFluidFieldList(Tensor::zero, "sigma normalization");
  FieldList<Dimension, Scalar> gdvNormalization = dataBase.newFluidFieldList(0.0, "grad div v normalization");

  // Grab the connectivity map from the DataBase.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const int numNodeLists = dataBase.numFluidNodeLists();
  CHECK(nodeLists.size() == numNodeLists);

  // Iterate over the NodeLists.
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // State for node i.
      const Vector& ri = position(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const SymTensor& Hi = Hfield(nodeListi, i);
      const Scalar weighti = mass(nodeListi, i)/rho(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // The derivatives for i.
      Tensor& sigmai = mSigma(nodeListi, i);
      Vector& gradDivVelocityi = mGradDivVelocity(nodeListi, i);
      Tensor& sigNormalizationi = sigNormalization(nodeListi, i);
      Scalar& gdvNormalizationi = gdvNormalization(nodeListi, i);

      // Connectivity for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Iterate over the neighboring NodeLists.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors in this NodeList.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            CHECK(j < nodeLists[nodeListj]->numNodes());

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {

              // State for node j.
              const Vector& rj = position(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const SymTensor& Hj = Hfield(nodeListj, j);
              const Scalar weightj = mass(nodeListj, j)/rho(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();

              // The derivatives for j.
              Tensor& sigmaj = mSigma(nodeListj, j);
              Vector& gradDivVelocityj = mGradDivVelocity(nodeListj, j);
              Tensor& sigNormalizationj = sigNormalization(nodeListj, j);
              Scalar& gdvNormalizationj = gdvNormalization(nodeListj, j);

              // Get an estimate of the smoothing scale for i and j.
              CHECK(distinctlyGreaterThan(Hdeti, 0.0));
              CHECK(distinctlyGreaterThan(Hdetj, 0.0));
              const Scalar frach = 0.02/(Dimension::rootnu(Hdeti) + Dimension::rootnu(Hdetj));
              const Scalar frach2 = frach*frach;
              CHECK(distinctlyGreaterThan(frach, 0.0));
              CHECK(distinctlyGreaterThan(frach2, 0.0));

              // Compute the relative distance and eta for i & j.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Scalar etamagi = etai.magnitude();
              const Vector Hetai = Hi*etai.unitVector();
              const pair<double, double> WWi = W.kernelAndGradValue(etamagi, Hdeti);
              const Scalar Wi = WWi.first;
              const Vector gWi = Hetai*WWi.second;

              const Vector etaj = Hj*rij;
              const Scalar etamagj = etaj.magnitude();
              const Vector Hetaj = Hj*etaj.unitVector();
              const pair<double, double> WWj = W.kernelAndGradValue(etamagj, Hdetj);
              const Scalar Wj = WWj.first;
              const Vector gWj = Hetaj*WWj.second;

              // Sum this pairs contribution to grad v for both i & j.
              const Vector vij = vi - vj;

              const double wij = weightj*Wi;
              const Vector gwij = weightj*gWi;

              const double wji = weighti*Wj;
              const Vector gwji = -weighti*gWj;

              const Vector rjiUnit = -rij.unitVector();
              const Tensor R = rotationMatrix(rjiUnit);
              const Tensor Rinverse = R.Transpose();
              const Scalar rij2 = rij.magnitude2();
              const Scalar rijmag = sqrt(rij2);

              const Scalar dxp = rijmag/(rij2 + frach2);
              const Vector dfdxpcol = -dxp*(R*vij);
              Tensor dfdxp;
              dfdxp.setColumn(0, dfdxpcol);
              dfdxp.rotationalTransform(Rinverse);

              sigmai += wij*dfdxp;
              sigmaj += wji*dfdxp;
            
              // Update the sigma normalization.
              Tensor tweighti = constructTensorWithColumnValue<Tensor, 0>(wij);
              Tensor tweightj = constructTensorWithColumnValue<Tensor, 0>(wji);
              tweighti.rotationalTransform(Rinverse);
              tweightj.rotationalTransform(-Rinverse);
              inPlaceAbsAdd<Tensor>(sigNormalizationi, tweighti);
              inPlaceAbsAdd<Tensor>(sigNormalizationj, tweightj);

              // Sum this pairs contribution to the second derivative of v.
              const Scalar thpt = vij.dot(rij)/(rij2 + frach2);

              gradDivVelocityi += thpt*gwij;
              gdvNormalizationi += wij;

              gradDivVelocityj += thpt*gwji;
              gdvNormalizationj += wji;

            }
          }
        }
      }

      // Complete the sigma calculation for i.
      CHECK(*std::min_element(sigNormalizationi.begin(), sigNormalizationi.end()) >= 0.0);
      tensorElementWiseDivide<Tensor>(sigmai, sigNormalizationi + tiny*Tensor::one);

      // Now limit to just negative eigen-values.  This is 'cause we only
      // care about convergent geometries for the Q.
      const SymTensor sigmai_s = sigmai.Symmetric();
      const Tensor sigmai_a = sigmai.SkewSymmetric();
      typename SymTensor::EigenStructType eigeni = sigmai_s.eigenVectors();
      sigmai = constructTensorWithMinDiagonal(eigeni.eigenValues, 0.0);
      sigmai.rotationalTransform(eigeni.eigenVectors);
      sigmai += sigmai_a;

      // Complete grad div v.
      CHECK(gdvNormalizationi >= 0.0);
      gradDivVelocityi /= gdvNormalizationi + tiny;

    }
  }

  // Apply boundary conditions
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
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
}
