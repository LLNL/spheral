//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructred the 
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#include "FiniteVolumeViscosity.hh"
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
#include "Mesh/Mesh.hh"
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
FiniteVolumeViscosity<Dimension>::
FiniteVolumeViscosity(const Scalar Clinear,
                      const Scalar Cquadratic,
                      const TableKernel<Dimension>& WT):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic, WT) {
}

//------------------------------------------------------------------------------
// Add our time derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
FiniteVolumeViscosity<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("FiniteVolumeViscosity_evalDerivs")

  // A few useful constants
  const auto Cl = this->mClinear;
  const auto Cq = this->mCquadratic;
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


      // Compute the pair QPi
      const auto vij = vi - vj;
      const auto xji = xj - xi;
      const auto xjihat = xji.unitVector();
      const auto hi = 1.0/(Hi*xjihat).magnitude();
      const auto hj = 1.0/(Hj*xjihat).magnitude();
      const auto DvDxi = min(0.0, DvDx(nodeListi, i).Trace());
      const auto DvDxj = min(0.0, DvDx(nodeListj, j).Trace());
      QPi[kk].first  = (-Clij*ci*DvDxi + Cqij*hi*DvDxi*DvDxi)*hi/rhoi;
      QPi[kk].second = (-Clij*cj*DvDxj + Cqij*hj*DvDxj*DvDxj)*hj/rhoj;

      // Stuff for time step constraints
      maxViscousPressurei = std::max(maxViscousPressurei, rhoi*rhoi*QPi[kk].first);
      maxViscousPressurej = std::max(maxViscousPressurej, rhoj*rhoj*QPi[kk].second);
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }      // OpenMP parallel region

  TIME_END("FiniteVolumeViscosity_evalDerivs")
}

//------------------------------------------------------------------------------
// Override the base method of computing the velocity gradient with a
// finite-volume approach
//------------------------------------------------------------------------------
template<typename Dimension>
void
FiniteVolumeViscosity<Dimension>::
updateVelocityGradient(const DataBase<Dimension>& dataBase,
                       const State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("FiniteVolumeViscosity_updateVelocityGradient")

  using Zone = typename Mesh<Dimension>::Zone;
  using Face = typename Mesh<Dimension>::Face;

  // Grab the DvDx for updating
  auto DvDx = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero);
  DvDx.Zero();

  // Make a finite-volume estimate of the local (to each Voronoi cell) velocity
  // gradient.
  unsigned nodeListj, j;
  Scalar Vi;
  const Mesh<Dimension>& mesh = state.mesh();
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const unsigned numNodeLists = velocity.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = velocity[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Vector& vi = velocity(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      const Zone& zonei = mesh.zone(nodeListi, i);
      const vector<int>& faces = zonei.faceIDs();
      Vi = zonei.volume();
      for (vector<int>::const_iterator fitr = faces.begin();
           fitr != faces.end();
           ++fitr) {
        const Face& faceij = mesh.face(*fitr);
        const int oppZoneID = faceij.oppositeZoneID(zonei.ID());
        if (Mesh<Dimension>::positiveID(oppZoneID) == Mesh<Dimension>::UNSETID) {
          nodeListj = nodeListi;
          j = i;
        } else {
          mesh.lookupNodeListID(Mesh<Dimension>::positiveID(oppZoneID), nodeListj, j);
        }
        const Vector& vj = velocity(nodeListj, j);
        const Vector vij = 0.5*(vi + vj);
        const Vector dA = faceij.area() * faceij.unitNormal() * sgn(*fitr);
        DvDxi -= vij*dA;
      }
      DvDxi /= Vi;
    }
  }
  TIME_END("FiniteVolumeViscosity_updateVelocityGradient")
}

}
