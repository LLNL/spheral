//---------------------------------Spheral++----------------------------------//
// ProbabilisticDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent scalar damage state.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//   Randles, P.W. & Libersky, L.D., 1996, "Comput. Methods Appl. Engrg, 
//     139, 375-408
//
// Created by JMO, Sun Sep 26 16:15:17 PDT 2004
//----------------------------------------------------------------------------//
#include "ProbabilisticDamagePolicy.hh"
#include "ProbabilisticDamageModel.hh"
#include "oneMinusDamage.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
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

namespace {

//------------------------------------------------------------------------------
// Force all eigen values positive.
//------------------------------------------------------------------------------
inline
void
abs_in_place(Dim<1>::Vector& vec) {
  vec[0] = std::abs(vec[0]);
}

inline
void
abs_in_place(Dim<2>::Vector& vec) {
  vec[0] = std::abs(vec[0]);
  vec[1] = std::abs(vec[1]);
}

inline
void
abs_in_place(Dim<3>::Vector& vec) {
  vec[0] = std::abs(vec[0]);
  vec[1] = std::abs(vec[1]);
  vec[2] = std::abs(vec[2]);
}

//------------------------------------------------------------------------------
// Sort the eigen values (& associated eigen vectors) in decreasing order.
//------------------------------------------------------------------------------
inline
void
sortEigen(Dim<1>::SymTensor::EigenStructType&) {
}

inline
void
sortEigen(Dim<2>::SymTensor::EigenStructType& eigeni) {
  if (eigeni.eigenValues(0) < eigeni.eigenValues(1)) {
    std::swap(eigeni.eigenValues(0), eigeni.eigenValues(1));
    eigeni.eigenVectors = Dim<2>::Tensor(eigeni.eigenVectors.xy(), eigeni.eigenVectors.xx(),
                                         eigeni.eigenVectors.yy(), eigeni.eigenVectors.yx());
  }
  ENSURE(eigeni.eigenValues(0) >= eigeni.eigenValues(1));
}

inline
void
sortEigen(Dim<3>::SymTensor::EigenStructType& eigeni) {
  int i = 0;
  int j = 1;
  int k = 2;
  bool flipped = true;
  while (flipped) {
    flipped = false;
    if (eigeni.eigenValues(i) < eigeni.eigenValues(j)) {
      flipped = true;
      std::swap(i, j);
    }
    if (eigeni.eigenValues(j) < eigeni.eigenValues(k)) {
      flipped = true;
      std::swap(j, k);
    }
  }
  eigeni.eigenValues = Dim<3>::Vector(eigeni.eigenValues(i),
                                      eigeni.eigenValues(j),
                                      eigeni.eigenValues(k));
  eigeni.eigenVectors = Dim<3>::Tensor(eigeni.eigenVectors(0,i), eigeni.eigenVectors(0,j), eigeni.eigenVectors(0,k),
                                       eigeni.eigenVectors(1,i), eigeni.eigenVectors(1,j), eigeni.eigenVectors(1,k),
                                       eigeni.eigenVectors(2,i), eigeni.eigenVectors(2,j), eigeni.eigenVectors(2,k));
  ENSURE((eigeni.eigenValues(0) >= eigeni.eigenValues(1)) &&
         (eigeni.eigenValues(1) >= eigeni.eigenValues(2)));
}

//------------------------------------------------------------------------------
// Determine the effective rotational transformation from the given velocity
// gradient.
//------------------------------------------------------------------------------
inline
Dim<1>::Tensor
effectiveRotation(const Dim<1>::Tensor&) {
  return Dim<1>::Tensor::one;
}

inline
Dim<2>::Tensor
effectiveRotation(const Dim<2>::Tensor& DvDx) {
  const double theta = DvDx.xy() - DvDx.yx();
  return Dim<2>::Tensor(cos(theta), sin(theta),
                        -sin(theta), cos(theta));
}

inline
Dim<3>::Tensor
effectiveRotation(const Dim<3>::Tensor& DvDx) {
  const Dim<3>::Vector theta(DvDx.yz() - DvDx.zy(),
                             DvDx.zx() - DvDx.xz(),
                             DvDx.xy() - DvDx.yx());
  return (Dim<3>::Tensor(cos(theta.x()), sin(theta.x()), 0.0,
                         -sin(theta.x()), cos(theta.x()), 0.0,
                         0.0, 0.0, 1.0)*
          Dim<3>::Tensor(cos(theta.y()), sin(theta.y()), 0.0,
                         0.0, 1.0, 0.0,
                         -sin(theta.y()), 0.0, cos(theta.y()))*
          Dim<3>::Tensor(1.0, 0.0, 0.0,
                         0.0, cos(theta.z()), sin(theta.z()),
                         0.0, -sin(theta.z()), cos(theta.z())));
}

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ProbabilisticDamagePolicy<Dimension>::
ProbabilisticDamagePolicy(const bool damageInCompression,  // allow damage in compression
                          const double kWeibull,           // coefficient in Weibull power-law
                          const double mWeibull,           // exponenent in Weibull power-law
                          const size_t minFlawsPerNode,    // minimum number of flaws to seed on any node
                          const double Vmin,               // minimum (initial) node volume
                          const double Vmax):              // maximum (initial) node volume
  UpdatePolicyBase<Dimension>(SolidFieldNames::strain),
  mDamageInCompression(damageInCompression),
  mMinFlawsPerNode(minFlawsPerNode),
  mkWeibull(kWeibull),
  mmWeibull(mWeibull),
  mVmin(Vmin),
  mVmax(Vmax) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ProbabilisticDamagePolicy<Dimension>::
~ProbabilisticDamagePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamagePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::tensorDamage);
  auto& stateField = state.field(key, SymTensor::zero);

  //const double tiny = 1.0e-30;
  //const double tol = 1.0e-5;

  // Get the state fields.
  const auto strainKey = State<Dimension>::buildFieldKey(SolidFieldNames::effectiveStrainTensor, nodeListKey);
  const auto DdamageDtKey = State<Dimension>::buildFieldKey(this->prefix() + SolidFieldNames::scalarDamage, nodeListKey);
  const auto DvDxKey = State<Dimension>::buildFieldKey(HydroFieldNames::internalVelocityGradient, nodeListKey);
  const auto numFlawsKey = State<Dimension>::buildFieldKey(SolidFieldNames::numFlaws, nodeListKey);
  const auto numFlawsActivatedKey = State<Dimension>::buildFieldKey(SolidFieldNames::numFlawsActivated, nodeListKey);
  const auto currentFlawKey = State<Dimension>::buildFieldKey(SolidFieldNames::currentFlaw, nodeListKey);
  const auto initialVolumeKey = State<Dimension>::buildFieldKey(SolidFieldNames::initialVolume, nodeListKey);
  const auto randomGeneratorsKey = State<Dimension>::buildFieldKey(SolidFieldNames::randomGenerators, nodeListKey);
  CHECK(state.registered(strainKey));
  CHECK(derivs.registered(DdamageDtKey));
  CHECK(derivs.registered(DvDxKey));
  CHECK(state.registered(numFlawsKey));
  CHECK(state.registered(numFlawsActivatedKey));
  CHECK(state.registered(currentFlawKey));
  CHECK(state.registered(initialVolumeKey));
  CHECK(state.registered(randomGeneratorsKey));
  const auto& strain = state.field(strainKey, SymTensor::zero);
  const auto& DDDt = derivs.field(DdamageDtKey, 0.0);
  const auto& localDvDx = derivs.field(DvDxKey, Tensor::zero);
  const auto& numFlaws = state.field(numFlawsKey, 0u);
  auto&       numFlawsActivated = state.field(numFlawsActivatedKey, 0u);
  auto&       currentFlaw = state.field(currentFlawKey, 0.0);
  const auto& initialVolume = state.field(initialVolumeKey, 0.0);
  auto&       randomGenerators = state.field(randomGeneratorsKey, uniform_random_01());

  // Iterate over the internal nodes.
  const auto ni = stateField.numInternalElements();
  const auto mInv = 1.0/mmWeibull;
#pragma omp parallel for
  for (auto i = 0u; i < ni; ++i) {
    auto& Di = stateField(i);
    CHECK(Di.eigenValues().minElement() >= 0.0 and
          fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));
    
    // Are we damaging this point?
    if (numFlaws(i) > 0u) {

      // First apply the rotational term to the current damage.
      const auto spin = localDvDx(i).SkewSymmetric();
      const auto spinCorrection = (spin*Di - Di*spin).Symmetric();
      Di += multiplier*spinCorrection;
      Di = max(1.0e-5, min(1.0 - 2.0e-5, Di));
      {
        const auto maxValue = Di.eigenValues().maxElement();
        if (maxValue > 1.0) Di /= maxValue;
      }
      CHECK(Di.eigenValues().minElement() >= 0.0 and
            fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));

      // The tensor strain on this node.
      auto eigeni = strain(i).eigenVectors();

      // If we're allowing damage in compression, force all strains to be abs value.
      if (mDamageInCompression) abs_in_place(eigeni.eigenValues);

      // We want to go over these from max to min.
      sortEigen(eigeni);

      // The direction of maximum strain absorbs the damage on this iteration.  This means we
      // only need the first of our sorted strain eigenvalues.
      // Only increment the damage in this direction if there is a positive strain.
      const auto straini = eigeni.eigenValues(0); //   + plasticStrain(i))/(fDi*fDi + 1.0e-20);

      // Are we activating new flaws?
      if (numFlawsActivated(i) < numFlaws(i) and straini > currentFlaw(i)) {
        ++numFlawsActivated(i);
        CHECK(numFlawsActivated(i) <= numFlaws(i));

        // The strain exceeds the current flaw activation threshold, so look for the next
        // value in the sequence, until we run out of flaws.
        if (numFlawsActivated(i) < numFlaws(i)) {
          vector<double> flaws(numFlaws(i) - numFlawsActivated(i));
          const auto Ai = numFlaws(i)/(mkWeibull*initialVolume(i));
          for (auto j = numFlawsActivated(i); j < numFlaws(i); ++j) {
            auto flawj = pow(Ai * randomGenerators(i)(), mInv);
            while (flawj < currentFlaw(i)) flawj = pow(Ai * randomGenerators(i)(), mInv);
            flaws[j] = flawj;
          }
          std::sort(flaws.begin(), flaws.end());

          // How many of the new flaws are activated?
          const auto itr = std::lower_bound(flaws.begin(), flaws.end(), straini);
          numFlawsActivated(i) += std::distance(flaws.begin(), itr);
          CHECK(numFlawsActivated(i) <= numFlaws(i));
          if (itr == flaws.end()) {
            currentFlaw(i) = flaws.back();
          } else {
            currentFlaw(i) = *itr;
          }
        }
      }

      // The direction of the strain, and projected current (starting) damage.
      const auto strainDirection = eigeni.eigenVectors.getColumn(0);
      CHECK(fuzzyEqual(strainDirection.magnitude2(), 1.0, 1.0e-10));
      const auto D0 = max(0.0, min(1.0, (Di * strainDirection).magnitude()));
      CHECK(D0 >= 0.0 && D0 <= 1.0);

      // Choose the allowed range of D.
      double Dmin, Dmax;
      const double numFlawsInv = 1.0/double(std::max(1u, numFlaws(i)));
      if (multiplier >= 0.0) {
        Dmin = D0;
        Dmax = std::max(D0, double(numFlawsActivated(i))*numFlawsInv);
      } else {
        Dmin = 0.0;
        Dmax = D0;
      }
      CHECK(Dmax >= Dmin);
      CHECK(Dmin >= 0.0 && Dmax <= 1.0);

      // Increment the damage.
      const auto D013 = FastMath::CubeRootHalley2(D0);
      const auto D113 = D013 + multiplier*double(numFlawsActivated(i))*numFlawsInv*DDDt(i);
      const auto D1 = max(Dmin, min(Dmax, FastMath::cube(D113)));
      CHECK((D1 - D0)*multiplier >= 0.0);

      // Increment the damage tensor.
      Di += (D1 - D0)*strainDirection.selfdyad();
    }

    // Enforce bounds on the damage.
    const auto Dvals = Di.eigenValues();
    if (Dvals.minElement() < 0.0 or
        Dvals.maxElement() > 1.0) {
      const auto Deigen = Di.eigenVectors();
      Di = constructSymTensorWithBoundedDiagonal(Deigen.eigenValues,
                                                 1.0e-5,
                                                 1.0);
      Di.rotationalTransform(Deigen.eigenVectors);
    }
    ENSURE(Di.eigenValues().minElement() >= 0.0 and
           fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));
    //     ENSURE(fuzzyGreaterThanOrEqual(Di.eigenValues().minElement(), 0.0, 1.0e-5) &&
    //            fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));
    //     ENSURE(fuzzyGreaterThanOrEqual(Di.eigenValues().minElement(), 0.0) &&
    //            fuzzyLessThanOrEqual(Di.Trace(), 1.0));

  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ProbabilisticDamagePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ProbabilisticDamagePolicy<Dimension>* rhsPtr = dynamic_cast<const ProbabilisticDamagePolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

