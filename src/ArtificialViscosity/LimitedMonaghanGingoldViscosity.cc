//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#include "LimitedMonaghanGingoldViscosity.hh"
#include "Boundary/Boundary.hh"
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

namespace {

//------------------------------------------------------------------------------
// limiter for velocity projection.
//------------------------------------------------------------------------------
//double limiterBJ(const double x) {
//  if (x > 0.0) {
//    return min(1.0, 4.0/(x + 1.0)*min(1.0, x));  // Barth-Jesperson
//  } else {
//    return 0.0;
//  }
//}
//
//double limiterMC(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*min(2.0*x, min(0.5*(1.0 + x), 2.0));   // monotonized central
//  } else {
//    return 0.0;
//  }
//}

double limiterVL(const double x) {
  if (x > 0.0) {
    return 2.0/(1.0 + x)*2.0*x/(1.0 + x);                       // van Leer
  } else {
    return 0.0;
  }
}

//double limiterMM(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*min(1.0, x);                           // minmod
//  } else {
//    return 0.0;
//  }
//}
//
//double limiterSB(const double x) {
//  if (x > 0.0) {
//    return 2.0/(1.0 + x)*max(min(2.0*x, 1.0), min(x, 2.0));    // superbee
//  } else {
//    return 0.0;
//  }
//}

// template<typename Vector, typename Tensor>
// double limiterConservative(const Vector& vi, const Vector& vj,
//                            const Vector& xi, const Vector& xj,
//                            const Tensor& DvDxi, const Tensor& DvDxj) {
//   const auto xji = xj - xi;
//   const auto vji = vj - vi;
//   const auto di = DvDxi.dot(xji);
//   const auto dj = DvDxj.dot(xji);
//   if (di.dot(dj) <= 0.0 or
//       di.dot(vji) <= 0.0 or
//       dj.dot(vji) <= 0.0) return 0.0;
//   const auto vjimag = vji.magnitude();
//   const auto dimag = di.magnitude();
//   const auto djmag = dj.magnitude();
//   return min(1.0, min(abs(vjimag*safeInv(dimag)), abs(vjimag*safeInv(djmag))));
// }

//double limiterConservative(const double vji, const double deltavi, const double deltavj) {
//  if (deltavi*deltavj <= 0.0 or
//      deltavi*vji <= 0.0 or
//      deltavj*vji <= 0.0) return 0.0;
//  return min(1.0, min(abs(vji*safeInv(deltavi)), abs(vji*safeInv(deltavj))));
//}

}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
LimitedMonaghanGingoldViscosity<Dimension>::
LimitedMonaghanGingoldViscosity(const Scalar Clinear,
                                const Scalar Cquadratic,
                                const bool linearInExpansion,
                                const bool quadraticInExpansion,
                                const Scalar etaCritFrac,
                                const Scalar etaFoldFrac,
                                const TableKernel<Dimension>& kernel):
  MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic, 
                                      linearInExpansion, quadraticInExpansion, kernel),
  mEtaCritFrac(etaCritFrac),
  mEtaFoldFrac(etaFoldFrac) {
}

//------------------------------------------------------------------------------
// Add our time derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
LimitedMonaghanGingoldViscosity<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("LimitedMonaghanGingoldViscosity_evalDerivs");

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

  // We need nPerh to figure out our critical folding distance. We assume the first NodeList value for this is
  // correct for all of them...
  const auto nPerh = position[0]->nodeList().nodesPerSmoothingScale();
  const auto etaCrit = mEtaCritFrac/nPerh;
  const auto etaFold = mEtaFoldFrac/nPerh;
  CHECK(etaFold > 0.0);

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
      const auto& DvDxi = DvDx(nodeListi, i);
      auto&       maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);

      const auto& xj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  fClj = noClCqMult ? 1.0 : ClMultiplier(nodeListj, j);
      const auto  fCqj = noClCqMult ? 1.0 : CqMultiplier(nodeListj, j);
      const auto& DvDxj = DvDx(nodeListj, j);
      auto&       maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);

      // Find the locally scaled coefficients
      const auto fshear = (balsaraCorrection ?
                           0.5*(this->calcBalsaraShearCorrection(DvDxi, Hi, ci) +
                                this->calcBalsaraShearCorrection(DvDxj, Hj, cj)) :
                           1.0);
      const auto Clij = 0.5*(fCli + fClj)*fshear * Cl;
      const auto Cqij = 0.5*(fCqi + fCqj)*fshear * Cq;

      // Displacement
      const auto xij = xi - xj;
      const auto etai = Hi*xij;
      const auto etaj = Hj*xij;

      // Compute the corrected velocity difference.
      const auto gradi = (DvDxi.dot(xij)).dot(xij);
      const auto gradj = (DvDxj.dot(xij)).dot(xij);
      const auto ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
      const auto rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
      CHECK(min(ri, rj) <= 1.0);
      // const Scalar phi = limiterMM(min(ri, rj));
      auto phi = limiterVL(min(ri, rj));
  
      // const auto xjihat = -xij.unitVector();
      // auto phi = limiterConservative((vj - vi).dot(xjihat), (DvDxi*xjihat).dot(xjihat), (DvDxj*xjihat).dot(xjihat));

      // If the points are getting too close, we let the Q come back full force.
      const auto etaij = min(etai.magnitude(), etaj.magnitude());
      // phi *= (etaij2 < etaCrit2 ? 0.0 : 1.0);
      // phi *= min(1.0, etaij2*etaij2/(etaCrit2etaCrit2));
      if (etaij < etaCrit) {
        phi *= exp(-FastMath::square((etaij - etaCrit)/etaFold));
      }

      // "Mike" method.
      const auto vi1 = vi - phi*DvDxi*xij;
      const auto vj1 = vj + phi*DvDxj*xij;

      // const Vector vi1 = vi - phi*DvDxi*xij;
      // const Vector vj1 = vj + phi*DvDxj*xij;
  
      const auto vij = vi1 - vj1;
  
      // Compute mu.
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

  TIME_END("LimitedMonaghanGingoldViscosity_evalDerivs");
}

}
