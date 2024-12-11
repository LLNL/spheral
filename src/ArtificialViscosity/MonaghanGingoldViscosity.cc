//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
//----------------------------------------------------------------------------//
#include "MonaghanGingoldViscosity.hh"
#include "DataOutput/Restart.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Neighbor/PairwiseField.hh"
#include "Material/EquationOfState.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "Utilities/Timer.hh"

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
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MonaghanGingoldViscosity<Dimension>::
MonaghanGingoldViscosity(const Scalar Clinear,
                         const Scalar Cquadratic,
                         const bool linearInExpansion,
                         const bool quadraticInExpansion,
                         const TableKernel<Dimension>& kernel):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic, kernel),
  mLinearInExpansion(linearInExpansion),
  mQuadraticInExpansion(quadraticInExpansion),
  mPairQPiPtr() {
}

//------------------------------------------------------------------------------
// Register derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
MonaghanGingoldViscosity<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  ArtificialViscosity<Dimension>::registerDerivatives(dataBase, derivs);
  const auto& connectivityMap = dataBase.connectivityMap();
  mPairQPiPtr = std::make_unique<PairQPiType>(connectivityMap);
  derivs.enroll(HydroFieldNames::pairQPi, *mPairQPiPtr);
}

//------------------------------------------------------------------------------
// Add our time derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
MonaghanGingoldViscosity<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("MonaghanGingoldViscosity_evalDerivs")

  // A few useful constants
  const auto Cl = this->mClinear;
  const auto Cq = this->mCquadratic;
  const auto eps2 = this->mEpsilon2;
  const auto balsaraCorrection = this->balsaraShearCorrection();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto ClMultiplier = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
  const auto CqMultiplier = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  const auto DvDx = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(ClMultiplier.size() == CqMultiplier.size());
  CHECK(DvDx.size() == numNodeLists);

  // Derivative FieldLists.
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto& QPi = derivs.template get<PairQPiType>(HydroFieldNames::pairQPi);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(QPi.size() == npairs);

  // Check if someone is evolving Cl and Cq coefficients
  const auto noClCqMult = ClMultiplier.size() == 0u;
  CHECK(noClCqMult or ClMultiplier.size() == numNodeLists);

  // Walk all the interacting pairs.
#pragma omp parallel
  {

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto i = pairs[kk].i_node;
      const auto j = pairs[kk].j_node;
      const auto nodeListi = pairs[kk].i_list;
      const auto nodeListj = pairs[kk].j_list;

      const auto& xi = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  fCli = noClCqMult ? 1.0 : ClMultiplier(nodeListi, i);
      const auto  fCqi = noClCqMult ? 1.0 : CqMultiplier(nodeListi, i);
      auto&       maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);

      const auto& xj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  fClj = noClCqMult ? 1.0 : ClMultiplier(nodeListj, j);
      const auto  fCqj = noClCqMult ? 1.0 : CqMultiplier(nodeListj, j);
      auto&       maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);

      // Find the locally scaled coefficients
      const auto fshear = (balsaraCorrection ?
                           0.5*(this->calcBalsaraShearCorrection(DvDx(nodeListi, i), Hi, ci) +
                                this->calcBalsaraShearCorrection(DvDx(nodeListj, j), Hj, cj)) :
                           1.0);
      const auto Clij = 0.5*(fCli + fClj)*fshear * Cl;
      const auto Cqij = 0.5*(fCqi + fCqj)*fshear * Cq;

      // Displacement
      const auto xij = xi - xj;
      const auto etai = Hi*xij;
      const auto etaj = Hj*xij;

      // Compute mu.
      const auto vij = vi - vj;
      const auto mui = vij.dot(etai)/(etai.magnitude2() + eps2);
      const auto muj = vij.dot(etaj)/(etaj.magnitude2() + eps2);

      // The artificial internal energy.
      const auto ei = -Clij*ci*(mLinearInExpansion    ? mui                : min(0.0, mui)) +
                       Cqij   *(mQuadraticInExpansion ? -sgn(mui)*mui*mui  : FastMath::square(min(0.0, mui)));
      const auto ej = -Clij*cj*(mLinearInExpansion    ? muj                : min(0.0, muj)) +
                       Cqij   *(mQuadraticInExpansion ? -sgn(muj)*muj*muj  : FastMath::square(min(0.0, muj)));
      CHECK2(ei >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ei << " " << ci << " " << mui);
      CHECK2(ej >= 0.0 or (mLinearInExpansion or mQuadraticInExpansion), ej << " " << cj << " " << muj);

      // Set the QPi value
      QPi[kk].first  = ei/rhoi;
      QPi[kk].second = ej/rhoj;

      // Stuff for time step constraints
      maxViscousPressurei = std::max(maxViscousPressurei, rhoi*rhoi*QPi[kk].first);
      maxViscousPressurej = std::max(maxViscousPressurej, rhoj*rhoj*QPi[kk].second);
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }      // OpenMP parallel region

  TIME_END("MonaghanGingoldViscosity_evalDerivs");
}

}
