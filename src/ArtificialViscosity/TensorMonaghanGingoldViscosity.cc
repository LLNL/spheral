//---------------------------------Spheral++----------------------------------//
// A modified form of the Monaghan & Gingold viscosity, extended to tensor 
// formalism.
//----------------------------------------------------------------------------//
#include "TensorMonaghanGingoldViscosity.hh"
#include "DataOutput/Restart.hh"
#include "Boundary/Boundary.hh"
#include "Geometry/EigenStruct.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Neighbor/PairwiseField.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/Timer.hh"
#include "Utilities/DBC.hh"

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
// Helper to remove any expansion terms from DvDx
//------------------------------------------------------------------------------
template<typename Tensor>
inline
void
removeExpansion(Tensor& DvDx) {
  const auto DvDx_s = DvDx.Symmetric();
  const auto DvDx_a = DvDx.SkewSymmetric();
  const auto eigeni = DvDx_s.eigenVectors();
  DvDx = constructTensorWithMinDiagonal(eigeni.eigenValues, 0.0);
  DvDx.rotationalTransform(eigeni.eigenVectors);
  DvDx += DvDx_a;
}

}

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorMonaghanGingoldViscosity<Dimension>::
TensorMonaghanGingoldViscosity(const Scalar Clinear,
                               const Scalar Cquadratic,
                               const TableKernel<Dimension>& kernel):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic, kernel),
  mPairQPiPtr() {
}

//------------------------------------------------------------------------------
// Register derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorMonaghanGingoldViscosity<Dimension>::
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
TensorMonaghanGingoldViscosity<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("TensorMonaghanGingoldViscosity_evalDerivs")

  // A few useful constants
  const auto tiny = 1.0e-20;
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
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto ClMultiplier = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);
  const auto CqMultiplier = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(ClMultiplier.size() == numNodeLists);
  CHECK(CqMultiplier.size() == numNodeLists);

  // Derivative FieldLists.
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto       maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto&      QPi = derivs.template get<PairQPiType>(HydroFieldNames::pairQPi);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(QPi.size() == npairs);

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Clij, Cqij, xij2, hi2, hj2, hi, hj, Qi, Qj, fshearij;
    Vector xij, vij, etai, etaj, xijUnit, thpt1, deltaSigmai, deltaSigmaj;
    Tensor R, Rinverse, mui, muj, Qepsi, Qepsj, QPiij, QPiji;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // If the nodes are not closing, then skip the rest and the Q for this pair is zero.
      const auto& xi = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& xj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      xij = xi - xj;
      vij = vi - vj;
      if (vij.dot(xij) < 0.0) {

        // Get the state for node i.
        const auto& rhoi = massDensity(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto& ci = soundSpeed(nodeListi, i);
        const auto& fCli = ClMultiplier(nodeListi, i);
        const auto& fCqi = CqMultiplier(nodeListi, i); 
        const auto& DvDxi = DvDx(nodeListi, i);
        auto&       maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
        CHECK(rhoi > 0.0);

        // Get the state for node j
        const auto& rhoj = massDensity(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto& cj = soundSpeed(nodeListj, j);
        const auto& fClj = ClMultiplier(nodeListj, j);
        const auto& fCqj = CqMultiplier(nodeListj, j); 
        const auto& DvDxj = DvDx(nodeListj, j);
        auto&       maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
        CHECK(rhoj > 0.0);

        // Node displacement.
        etai = Hi*xij;
        etaj = Hj*xij;

        // Find the effective linear and quadratic coefficients for this pair.
        // These incorporate things like the Morris & Monaghan time evolved coefficients.
        Clij = 0.5*(fCli + fClj) * Cl;
        Cqij = 0.5*(fCqi + fCqj) * Cq;

        // Scale by the Balsara shear correction if necessary
        if (balsaraCorrection) {
          fshearij = std::max(this->calcBalsaraShearCorrection(DvDxi, Hi, ci),
                              this->calcBalsaraShearCorrection(DvDxj, Hj, cj));
          Clij *= fshearij;
          Cqij *= fshearij;
        }

        // Some more geometry.
        xij2 = xij.magnitude2();
        xijUnit = xij.unitVector();
        hi2 = xij2/(etai.magnitude2() + tiny);
        hj2 = xij2/(etaj.magnitude2() + tiny);
        hi = sqrt(hi2);
        hj = sqrt(hj2);

        // Build the tensor for grad-v with the pair-wise value sliced in
        Tensor sigmai = DvDxi;
        Tensor sigmaj = DvDxj;
        {
          R = rotationMatrix(xijUnit);
          Rinverse = R.Transpose();
          thpt1 = sqrt(xij2)*(R*vij);
          deltaSigmai = thpt1/(xij2 + eps2*hi2);
          deltaSigmaj = thpt1/(xij2 + eps2*hj2);
          sigmai.rotationalTransform(R);
          sigmaj.rotationalTransform(R);
          sigmai.setColumn(0, deltaSigmai);
          sigmaj.setColumn(0, deltaSigmaj);
          sigmai.rotationalTransform(Rinverse);
          sigmaj.rotationalTransform(Rinverse);
        }

        // Remove any expansive components from sigma
        removeExpansion(sigmai);
        removeExpansion(sigmaj);

        // Calculate the tensor viscous internal energy.
        mui = hi*sigmai;
        Qepsi = -Clij*ci*mui.Transpose() + Cqij*mui*mui;

        muj = hj*sigmaj;
        Qepsj = -Clij*cj*muj.Transpose() + Cqij*muj*muj;

        // We now have enough to compute Pi!
        QPiij = Qepsi/rhoi;
        QPiji = Qepsj/rhoj;
        QPi[kk] = std::make_pair(QPiij, QPiji);

        // Stuff for time step constraints
        Qi = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
        Qj = rhoj*rhoj*(QPiji.diagonalElements().maxAbsElement());
        maxViscousPressurei = std::max(maxViscousPressurei, Qi);
        maxViscousPressurej = std::max(maxViscousPressurej, Qj);
      }  // If vij.dot(xij) < 0
    }    // pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }      // OpenMP parallel region

  TIME_END("TensorMonaghanGingoldViscosity_evalDerivs");
}

}

