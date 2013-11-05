//---------------------------------Spheral++----------------------------------//
// TensorDamageModel -- Base class for the tensor damage physics models.
// This class does not know how to seed the flaw distribution -- that is 
// required of descendant classes.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//   Randles, P.W. & Libersky, L.D., 1996, "Comput. Methods Appl. Engrg, 
//     139, 375-408
//
// Created by JMO, Thu Sep 29 15:42:05 PDT 2005
//----------------------------------------------------------------------------//
#include <string>
#include <vector>
#include <algorithm>

#include "boost/shared_ptr.hpp"

#include "TensorDamageModel.hh"
#include "TensorStrainPolicy.hh"
#include "TensorDamagePolicy.hh"
#include "EffectiveTensorDamagePolicy.hh"
#include "DamageGradientPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "Strength/SolidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/ReplaceState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FileIO/FileIO.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/Neighbor.hh"

namespace Spheral {
namespace PhysicsSpace {

using namespace std;

using NodeSpace::NodeList;
using SolidMaterial::SolidNodeList;
using Material::EquationOfState;
using FileIOSpace::FileIO;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using BoundarySpace::Boundary;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NeighborSpace::Neighbor;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorDamageModel<Dimension>::
TensorDamageModel(SolidNodeList<Dimension>& nodeList,
                  const TensorStrainAlgorithm strainAlgorithm,
                  const EffectiveDamageAlgorithm effDamageAlgorithm,
                  const bool useDamageGradient,
                  const TableKernel<Dimension>& W,
                  const double crackGrowthMultiplier,
                  const EffectiveFlawAlgorithm flawAlgorithm,
                  const double criticalDamageThreshold,
                  const FlawStorageType& flaws):
  DamageModel<Dimension>(nodeList, W, crackGrowthMultiplier, flawAlgorithm, flaws),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mEffectiveStrain(SolidFieldNames::effectiveStrainTensor, nodeList),
  mDdamageDt(TensorDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mNewEffectiveDamage(ReplaceState<Dimension, FieldSpace::Field<Dimension, SymTensor> >::prefix() + SolidFieldNames::effectiveTensorDamage, nodeList),
  mNewDamageGradient(ReplaceState<Dimension, FieldSpace::Field<Dimension, Vector> >::prefix() + SolidFieldNames::damageGradient, nodeList),
  mStrainAlgorithm(strainAlgorithm),
  mEffDamageAlgorithm(effDamageAlgorithm),
  mCriticalDamageThreshold(criticalDamageThreshold),
  mUseDamageGradient(useDamageGradient) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorDamageModel<Dimension>::
~TensorDamageModel() {
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
// For this model we evaluate the derivative of the tensor damage field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  const double tiny = 1.0e-15;

  // The base class determines the scalar magnitude of the damage evolution.
  typedef typename State<Dimension>::KeyType Key;
  const SolidNodeList<Dimension>* nodeListPtr = &(this->nodeList());
  Field<Dimension, Scalar>& DDDt = derivs.field(state.buildFieldKey(TensorDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeListPtr->name()), 0.0);
  this->computeScalarDDDt(dataBase,
                          state,
                          time,
                          dt,
                          DDDt);

  // We need to fill in the effective damage and perhaps the damage gradient.
  // First, grab the pertinent state.
  const Field<Dimension, Scalar>& mass = state.field(state.buildFieldKey(HydroFieldNames::mass, nodeListPtr->name()), 0.0);
  const Field<Dimension, Vector>& position = state.field(state.buildFieldKey(HydroFieldNames::position, nodeListPtr->name()), Vector::zero);
  const Field<Dimension, Scalar>& rho = state.field(state.buildFieldKey(HydroFieldNames::massDensity, nodeListPtr->name()), 0.0);
  const Field<Dimension, SymTensor>& H = state.field(state.buildFieldKey(HydroFieldNames::H, nodeListPtr->name()), SymTensor::zero);
  const Field<Dimension, SymTensor>& D = state.field(state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListPtr->name()), SymTensor::zero);

  Field<Dimension, SymTensor>& Deff = derivs.field(state.buildFieldKey(EffectiveTensorDamagePolicy<Dimension>::prefix() + SolidFieldNames::effectiveTensorDamage, nodeListPtr->name()), SymTensor::zero);
  Field<Dimension, Vector>& gradD = derivs.field(state.buildFieldKey(DamageGradientPolicy<Dimension>::prefix() + SolidFieldNames::damageGradient, nodeListPtr->name()), Vector::zero);

  // Get the kernel.
  const TableKernel<Dimension>& W = this->kernel();
  const double etaMax2 = W.kernelExtent() * W.kernelExtent();

  // The neighbor connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), nodeListPtr));
  CHECK(nodeListi < nodeLists.size());
  const int firstGhostNodei = nodeListPtr->firstGhostNode();

  // Prepare the effective damage for computation.
  Field<Dimension, Scalar> normalization("normalization", *nodeListPtr, 0.0);
  if (mEffDamageAlgorithm == Max or
      mEffDamageAlgorithm == Copy) Deff = D;

  // Iterate over the internal nodes in this NodeList.
  for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
       iItr != connectivityMap.end(nodeListi);
       ++iItr) {
    const int i = *iItr;

    // State for node i.
    const Scalar& mi = mass(i);
    const Vector& ri = position(i);
    const Scalar& rhoi = rho(i);
    const SymTensor& Hi = H(i);
    const Scalar Hdeti = Hi.Determinant();
    const SymTensor& Di = D(i);
    const Scalar Dmagi = max(0.0, min(1.0, Di.Trace()/Dimension::nDim));
    Scalar& normalizationi = normalization(i);
    SymTensor& Deffi = Deff(i);
    Vector& gradDi = gradD(i);
    CHECK(mi > 0.0);
    CHECK(rhoi > 0.0);
    CHECK(Hdeti > 0.0);
    CHECK(Dmagi >= 0.0 and Dmagi <= 1.0);

    // Get the connectivity info for this node.  We only need to proceed if
    // there are some nodes in this list.
    const vector<int>& connectivity = connectivityMap.connectivityForNode(nodeListPtr, i)[nodeListi];
    if (connectivity.size() > 0) {

      // Iterate over the neighbors.
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
        const int j = *jItr;

        // Only proceed if this node pair has not been calculated yet.
        if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                     nodeListi, j,
                                                     firstGhostNodei)) {

          // Get the state for node j
          const Scalar& mj = mass(j);
          const Vector& rj = position(j);
          const Scalar& rhoj = rho(j);
          const SymTensor& Hj = H(j);
          const Scalar Hdetj = Hj.Determinant();
          const SymTensor& Dj = D(j);
          const Scalar Dmagj = max(0.0, min(1.0, Dj.Trace()/Dimension::nDim));
          Scalar& normalizationj = normalization(j);
          SymTensor& Deffj = Deff(j);
          Vector& gradDj = gradD(j);
          CHECK(mj > 0.0);
          CHECK(rhoj > 0.0);
          CHECK(Hdetj > 0.0);
          CHECK(Dmagj >= 0.0 and Dmagj <= 1.0);

          // Node displacement and weighting.
          const Vector rij = ri - rj;

          const Vector etai = Hi*rij;
          const Scalar etaMagi = etai.magnitude();
          const std::pair<double, double> WWi = W.kernelAndGradValue(etaMagi, Hdeti);
          const Scalar Wi = WWi.first;
          const Vector gradWi = WWi.second * (Hi*etai.unitVector());

          const Vector etaj = Hj*rij;
          const Scalar etaMagj = etaj.magnitude();
          const std::pair<double, double> WWj = W.kernelAndGradValue(etaMagj, Hdetj);
          const Scalar Wj = WWj.first;
          const Vector gradWj = WWj.second * (Hj*etaj.unitVector());

          // Increment the effective damage.
          switch(mEffDamageAlgorithm) {
          case Max:
            if (Dmagj * Dimension::nDim > Deffi.Trace()) Deffi = Dj;
            if (Dmagi * Dimension::nDim > Deffj.Trace()) Deffj = Di;
            break;

          case Sampled:
            normalizationi += Wi;
            normalizationj += Wj;
            Deffi += Wi * Dj;
            Deffj += Wj * Di;
            break;
          }

          // Increment the gradient.
          gradDi += mj*(Dj - Di)*gradWi;
          gradDj += mi*(Dj - Di)*gradWj;

        }
      }
    }

    // Finish the effective damage by adding in the self-contribution and 
    // normalizing the sucker.
    if (mEffDamageAlgorithm == Sampled) {
      CHECK(normalizationi >= 0.0);
      const double Wi = W.kernelValue(0.0, Hdeti);
      Deffi = (Deffi + Dmagi*Wi*Di)/(normalizationi + Dmagi*Wi + tiny);
    }

    // Finish the gradient.
    if (mUseDamageGradient) {
      CHECK(rhoi > 0.0);
      gradDi /= rhoi;
    } else {
      gradDi.Zero();
    }

  }

}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename TensorDamageModel<Dimension>::TimeStepType
TensorDamageModel<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {

//   // Look at how quickly we're trying to change the damage.
//   double dt = DBL_MAX;
//   const Field<Dimension, SymTensor>& damage = this->nodeList().damage();
//   const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
//   const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
//   const size_t nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), &(this->nodeList())));
//   for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
//        iItr != connectivityMap.end(nodeListi);
//        ++iItr) {
//     const int i = *iItr;
//     const double D0 = damage(i).Trace() / Dimension::nDim;
//     dt = min(dt, 0.8*max(D0, 1.0 - D0)/
//              std::sqrt(mDdamageDt(i)*mDdamageDt(i) + 1.0e-20));
//   }

  return TimeStepType(1.0e100, "Rate of damage change -- NO VOTE.");
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::KeyType Key;
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Register the strain and effective strain.
  PolicyPointer effectiveStrainPolicy(new TensorStrainPolicy<Dimension>(mStrainAlgorithm));
  state.enroll(mStrain);
  state.enroll(mEffectiveStrain, effectiveStrainPolicy);

  // Register the damage and effective damage.
  // Note we are overriding the default no-op policy for the effective damage
  // as originally registered by the SolidSPHHydroBase class.
  PolicyPointer damagePolicy(new TensorDamagePolicy<Dimension>(*this));
  PolicyPointer effDamagePolicy(new EffectiveTensorDamagePolicy<Dimension>());
  state.enroll(this->nodeList().damage(), damagePolicy);
  state.enroll(this->nodeList().effectiveDamage(), effDamagePolicy);

  // If we're using the grad damage switch, override the policy for the updating
  // the damage gradient.
  if (mUseDamageGradient) {
    PolicyPointer damageGradientPolicy(new DamageGradientPolicy<Dimension>());
    state.enroll(this->nodeList().damageGradient(), damageGradientPolicy);
  }

  // // Mask out nodes beyond the critical damage threshold from setting the timestep.
  // Key maskKey = state.buildFieldKey(HydroFieldNames::timeStepMask, this->nodeList().name());
  // Field<Dimension, int>& mask = state.field(maskKey, 0);
  // const Field<Dimension, SymTensor>& damage = this->nodeList().damage();
  // for (int i = 0; i != this->nodeList().numInternalNodes(); ++i) {
  //   if (damage(i).Trace() > mCriticalDamageThreshold) mask(i) = 0;
  // }

//   // Damage some of the state variables.
//   typedef typename State<Dimension>::FieldKeyType Key;
//   const NodeList<Dimension>* nodeListPtr = &(this->nodeList());
//   const Key KKey(nodeListPtr, SolidFieldNames::bulkModulus);
//   const Key muKey(nodeListPtr, SolidFieldNames::shearModulus);
//   const Key DKey(nodeListPtr, SolidFieldNames::tensorDamage);
//   Field<Dimension, Scalar>& K = state.scalarField(KKey);
//   Field<Dimension, Scalar>& mu = state.scalarField(muKey);
//   const Field<Dimension, SymTensor>& D = state.symTensorField(DKey);
//   for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
//     const Scalar fDi = std::max(0.0, std::min(1.0, 1.0 - D(i).eigenValues().maxElement() - 1.0e-5));
//     CHECK(fDi >= 0.0 && fDi <= 1.0);
//     K(i) *= fDi;
//     mu(i) *= fDi;
//   }

  // Register the base classes stuff.
  DamageModel<Dimension>::registerState(dataBase, state);
}

//------------------------------------------------------------------------------
// Register the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDdamageDt);
  derivs.enroll(mNewEffectiveDamage);
  derivs.enroll(mNewDamageGradient);
}

//------------------------------------------------------------------------------
// Apply the boundary conditions to the ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Grab this models damage field from the state.
  typedef typename State<Dimension>::KeyType Key;
  const Key nodeListName = this->nodeList().name();
  const Key DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  const Key DeffKey = state.buildFieldKey(SolidFieldNames::effectiveTensorDamage, nodeListName);
  const Key gradDKey = state.buildFieldKey(SolidFieldNames::damageGradient, nodeListName);
  CHECK(state.registered(DKey));
  CHECK(state.registered(DeffKey));
  CHECK(state.registered(gradDKey));
  Field<Dimension, SymTensor>& D = state.field(DKey, SymTensor::zero);
  Field<Dimension, SymTensor>& Deff = state.field(DeffKey, SymTensor::zero);
  Field<Dimension, Vector>& gradD = state.field(gradDKey, Vector::zero);

  // Apply ghost boundaries to the damage.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyGhostBoundary(D);
    (*boundaryItr)->applyGhostBoundary(Deff);
    if (mUseDamageGradient) (*boundaryItr)->applyGhostBoundary(gradD);
  }
}

//------------------------------------------------------------------------------
// Enforce boundary conditions for the physics specific fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Grab this models damage field from the state.
  typedef typename State<Dimension>::KeyType Key;
  const Key nodeListName = this->nodeList().name();
  const Key DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  const Key DeffKey = state.buildFieldKey(SolidFieldNames::effectiveTensorDamage, nodeListName);
  const Key gradDKey = state.buildFieldKey(SolidFieldNames::damageGradient, nodeListName);
  CHECK(state.registered(DKey));
  CHECK(state.registered(DeffKey));
  CHECK(state.registered(gradDKey));
  Field<Dimension, SymTensor>& D = state.field(DKey, SymTensor::zero);
  Field<Dimension, SymTensor>& Deff = state.field(DeffKey, SymTensor::zero);
  Field<Dimension, Vector>& gradD = state.field(gradDKey, Vector::zero);

  // Enforce!
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceBoundary(D);
    (*boundaryItr)->enforceBoundary(Deff);
    if (mUseDamageGradient) (*boundaryItr)->enforceBoundary(gradD);
  }
}


// //------------------------------------------------------------------------------
// // Initialize at beginning of a step.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// TensorDamageModel<Dimension>::
// initialize(const typename Dimension::Scalar& time,
//            const typename Dimension::Scalar& dt,
//            const DataBase<Dimension>& dataBase,
//            State<Dimension>& state,
//            StateDerivatives<Dimension>& derivs) {

//   // We need to fill in the effective damage and perhaps the damage gradient.
//   // First, grab the pertinent state.
//   typedef typename State<Dimension>::FieldKeyType KeyType;
//   const SolidNodeList<Dimension>& nodeList = this->nodeList();
//   const Field<Dimension, Scalar>& mass = state.scalarField(KeyType(nodeListPtr, HydroFieldNames::mass));
//   const Field<Dimension, Vector>& position = state.vectorField(KeyType(nodeListPtr, HydroFieldNames::position));
//   const Field<Dimension, Scalar>& rho = state.scalarField(KeyType(nodeListPtr, HydroFieldNames::massDensity));
//   const Field<Dimension, Scalar>& weight = state.scalarField(KeyType(nodeListPtr, HydroFieldNames::weight));
//   const Field<Dimension, SymTensor>& H = state.symTensorField(KeyType(nodeListPtr, HydroFieldNames::H));
//   const Field<Dimension, Scalar>& omega = state.scalarField(KeyType(nodeListPtr, HydroFieldNames::omegaGradh));
//   const Field<Dimension, SymTensor>& D = state.symTensorField(KeyType(nodeListPtr, SolidFieldNames::tensorDamage));

//   Field<Dimension, SymTensor>& Deff = state.symTensorField(KeyType(nodeListPtr, SolidFieldNames::effectiveTensorDamage));
//   Field<Dimension, Vector>& gradD = derivs.vectorField(KeyType(&nodeList, SolidFieldNames::gradDamage));

//   // Get the kernel.
//   const TableKernel<Dimension>& W = nodeListPtr->kernel();
//   const double etaMax2 = W.kernelExtent() * W.kernelExtent();

//   // The neighbor connectivity.
//   const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
//   const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
//   const size_t nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), nodeListPtr));
//   CHECK(nodeListi < nodeLists.size());

//   const size_t firstGhostNode = nodeListPtr->firstGhostNode();

//   // Iterate over the internal nodes in this NodeList.
//   Field<Dimension, Scalar> normaliation("normalization", *nodeListPtr, 0.0);
//   Deff = D;
//   for (size_t i = 0; i != nodeListPtr->numInternalNodes(); ++i) {

//     // State for node i.
//     const Scalar& mi = mass(i);
//     const Vector& ri = position(i);
//     const Scalar& rhoi = rho(i);
//     const Scalar& weighti = weight(i);
//     const SymTensor& Hi = H(i);
//     const Scalar Hdeti = Hi.Determinant();
//     const Scalar safeOmegai = omega(i)/(omega(i)*omega(i) + 1.0e-20);
//     const SymTensor& Di = D(i);
//     const Scalar Dmagi = Di.Trace();
//     CHECK(mi > 0.0);
//     CHECK(rhoi > 0.0);
//     CHECK(weighti > 0.0);
//     CHECK(safeOmegai > 0.0);
//     CHECK(Hdeti > 0.0);
//     CHECK(Dmagi >= 0.0 and Dmagi <= Dimension::nDim);

//     SymTensor& Deffi = Deff(i);
//     Vector& gradDi = gradD(i);

//     // Get the connectivity info for this node.  We only need to proceed if
//     // there are some nodes in this list.
//     const vector<int>& connectivity = connectivityMap.connectivityForNode(nodeListPtr, i)[nodeListi];
//     if (connectivity.size() > 0) {

//       // Iterate over the neighbors.
//       for (vector<int>::const_iterator jItr = connectivity.begin();
//            jItr != connectivity.end();
//            ++jItr) {
//         const int j = *jItr;

//         // Only proceed if this node pair has not been calculated yet.
//         if (j > i) {

//           // Get the state for node j
//           const Scalar& mj = mass(j);
//           const Vector& rj = position(j);
//           const Scalar& rhoj = rho(j);
//           const Scalar& weightj = weight(j);
//           const SymTensor& Hj = H(j);
//           const Scalar Hdetj = Hj.Determinant();
//           const Scalar safeOmegaj = omega(j)/(omega(j)*omega(j) + 1.0e-20);
//           const SymTensor& Dj = D(j);
//           const Scalar& Dmagj = Dj.Trace();
//           CHECK(mj > 0.0);
//           CHECK(rhoj > 0.0);
//           CHECK(weightj > 0.0);
//           CHECK(safeOmegaj > 0.0);
//           CHECK(Hdetj > 0.0);
//           CHECK(Dmagj >= 0.0 and Dmagj <= Dimension::nDim);

//           SymTensor& Deffj = Deff(j);
//           Vector& gradDnewj = gradD(j);

//           // Node displacement and weighting.
//           const Vector rij = ri - rj;

//           const Vector etai = Hi*rij;
//           const Scalar etaMagi = etai.magnitude();
//           const std::pair<double, double> WWi = W.kernelAndGradValue(etaMagi, Hdeti);
//           const Scalar Wi = WWi.first;
//           const Vector gradWi = WWi.second * (Hi*etai.unitVector());

//           const Vector etaj = Hj*rij;
//           const Scalar etaMagj = etaj.magnitude();
//           const std::pair<double, double> WWj = W.kernelAndGradValue(etaMagj, Hdetj);
//           const Scalar Wj = WWj.first;
//           const Vector gradWj = WWj.second * (Hj*etaj.unitVector());


//           const Vector gradWj = W.grad(etaMagj, Hdetj) * (Hj*etaj.unitVector());

//           // Increment the effective damage.
//           switch(mEffDamageAlgorithm) {
//           case Max:
//             if (Dmagj > Dmaji) Di = Dj;
//             if (Dmagi > Dmajj) Dj = Di;
//             break;

//           case Sampled:
//             normalizationi += Dmagj*Wi;
//             normalizationj += Dmagi*Wj;
//             Di += Dmagj*Wi * Dj;
//             Dj += Dmagi*Wj * Di;
//             break;
//           }

//           // Increment the gradient.
//           gradDi += mj*(Dj - Di)*gradWi;
//           gradDj += mi*(Dj - Di)*gradWj;

//         }
//       }
//     }

//     // Finish the effective damage.
//     if (mEffDamageAlgorithm == Sampled) {
//       CHECK(normalizationi >= 0.0);
//       Deffi /= normalizationi + tiny;
//     }

//     // Finish the gradient.
//     if (mUseDamageGradient) {
//       CHECK(rhoi > 0.0);
//       gradDi /= rhoi;
//     } else {
//       gradDi.Zero();
//     }

//   }










//   // If needed, compute the gradient of the damage.
//   if (mUseDamageGradient) {

//     // State fields.
//     typedef typename State<Dimension>::FieldKeyType KeyType;
//     const SolidNodeList<Dimension>& nodeList = this->nodeList();
//     const Field<Dimension, Scalar>& mass = state.scalarField(KeyType(&nodeList, HydroFieldNames::mass));
//     const Field<Dimension, Vector>& position = state.vectorField(KeyType(&nodeList, HydroFieldNames::position));
//     const Field<Dimension, Scalar>& rho = state.scalarField(KeyType(&nodeList, HydroFieldNames::massDensity));
//     const Field<Dimension, SymTensor>& H = state.symTensorField(KeyType(&nodeList, HydroFieldNames::H));
//     const Field<Dimension, Scalar>& omega = state.scalarField(KeyType(nodeListPtr, HydroFieldNames::omegaGradh));
//     const Field<Dimension, Scalar>& D = state.scalarField(KeyType(&nodeList, SolidFieldNames::tensorDamage));
//     Field<Dimension, Vector>& gradD = derivs.vectorField(KeyType(&nodeList, SolidFieldNames::gradDamage));

//     // Get the kernel.
//     const TableKernel<Dimension>& W = nodeList.kernel();
//     const double etaMax2 = W.kernelExtent() * W.kernelExtent();

//     // The neighbor connectivity.
//     const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
//     const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
//     const int nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), &nodeList));
//     CHECK(nodeListi < nodeLists.size());

//     // Initialize the grad field.
//     gradD.Zero();

//     // Iterate over the nodes in this NodeList.
//     for (size_t i = 0; i != nodeList.numInternalNodes(); ++i) {

//       // State for node i.
//       const Vector& ri = position(i);
//       const Scalar& rhoi = rho(i);
//       const SymTensor& Hi = H(i);
//       const Scalar& Di = D(i);

//       // Get the connectivity for this node.
//       const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);
//       CHECK(fullConnectivity.size() == nodeLists.size());
//       const vector<int>& connectivity = fullConnectivity[nodeListi];

//       // Iterate over the neighbor nodes in this NodeList.
//       for (typename vector<int>::const_iterator jItr = connectivity.begin();
//            jItr != connectivity.end();
//            ++jItr) {
//         const int j = *jItr;
     
//         // State for node j.
//         const Scalar& mj = mass(j);
//         const Vector& rj = position(j);
//         const Scalar& Dj = D(j);

//         // Node displacement and weighting.
//         const Vector rij = ri - rj;
//         const Vector etai = Hi*rij;
//         const Vector gradWi = Hi*W.grad(etai.magnitude(), Hi.Determinant()) * etai.unitVector();

//         // Increment the gradient.
//         gradD(i) += mj*(Dj - Di)*gradWi;
//       }

//       // Finalize the gradient.
//       CHECK(rhoi > 0.0);
//       gradD(i) /= rhoi;

//     }

//     // Apply boundary conditions to the gradient.
//     for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
//          boundaryItr != this->boundaryEnd();
//          ++boundaryItr) {
//       (*boundaryItr)->applyGhostBoundary(gradD);
//     }
//     for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
//          boundaryItr != this->boundaryEnd();
//          ++boundaryItr) {
//       (*boundaryItr)->finalizeGhostBoundary();
//     }

//   }

// }

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  DamageModel<Dimension>::dumpState(file, pathName);
  file.write(mStrain, pathName + "/strain");
  file.write(mEffectiveStrain, pathName + "/effectiveStrain");
  file.write(mDdamageDt, pathName + "/DdamageDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamageModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  DamageModel<Dimension>::restoreState(file, pathName);
  file.read(mStrain, pathName + "/strain");
  file.read(mEffectiveStrain, pathName + "/effectiveStrain");
  file.read(mDdamageDt, pathName + "/DdamageDt");
}

}
}

