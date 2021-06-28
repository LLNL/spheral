//---------------------------------Spheral++----------------------------------//
// Physics -- dampled linear spring contact model Spheral++
//----------------------------------------------------------------------------//
#include "DampedLinearSpring.hh"

#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Hydro/HydroFieldNames.hh"

#include <cmath>

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
DampedLinearSpring<Dimension>::
DampedLinearSpring(Scalar YoungsModulus,
                   Scalar restitutionCoefficient):
                   mRestitutionCoefficient(restitutionCoefficient),
                   mYoungsModulus(YoungsModulus) {
      mBeta = 3.14159/std::log(restitutionCoefficient);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DampedLinearSpring<Dimension>::
~DampedLinearSpring() {}



template<typename Dimension>
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
timeStep(const DataBase<Dimension>& dataBase,
         const State<Dimension>& state,
         const StateDerivatives<Dimension>& /*derivs*/,
               typename Dimension::Scalar /*currentTime*/) const{
  const auto& mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const auto& mass = state.fields(HydroFieldNames::mass, 0.0); 
  const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto& velocity = state.fields(HydroFieldNames::velocity, Vector::zero);  
  const auto& angularVelocity = state.fields("angularVelocity", Vector::zero);  
  const auto& radius = state.fields("particleRadius",0.0);

  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  auto minContactTime = 1e30;
  const auto pi = 3.1415;
  #pragma omp parallel
  {

    int i, j, nodeListi, nodeListj;

  #pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {

      const auto mi = mass(nodeListi,i);
      const auto mj = mass(nodeListj,j);
      const auto Ri = radius(nodeListi,i);
      const auto Rj = radius(nodeListj,j);

      const auto mij = (mi*mj)/(mi+mj);
      const auto Rij = (Ri+Rj)/(Ri+Rj);

      const auto c1 = 4.0/3.0*mYoungsModulus*std::sqrt(Rij);
      const auto c2 = std::sqrt(4.0*mij*c1/(1+mBeta*mBeta));

      const auto stiffness = c1/mij;
      const auto dissipation = c2/(2.0*mij);
      const auto contactFrequency = std::sqrt(stiffness - dissipation*dissipation);
      const auto contactTime = pi/contactFrequency;
      minContactTime = min(contactTime,minContactTime);
    }
  }
  return minContactTime;
};

//------------------------------------------------------------------------------
// get our acceleration and other things
//------------------------------------------------------------------------------
template<typename Dimension>
void
DampedLinearSpring<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const{

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  
  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  //const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  //const auto omega = state.fields(Hydro, Vector::zero);
  const auto radius = state.fields("particleRadius",0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  //CHECK(H.size() == numNodeLists);
  //CHECK(omega.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DvDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);

  // The set of interacting node pairs.

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& Ri = radius(nodeListi, i);
      
      auto& DvDti = DvDt_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& Rj = radius(nodeListj, j);

      auto& DvDtj = DvDt_thread(nodeListj, j);

      CHECK(mi > 0.0);
      CHECK(Hdeti > 0.0);
      CHECK(mj > 0.0);
      CHECK(Hdetj > 0.0);

      const auto mij = (mi*mj)/(mi+mj);
      const auto Rij = (Ri*Rj)/(Ri+Rj);

      const auto vij = vi-vj;
      const auto rij = ri-rj;
      const auto rhatij = rij.unitVector();

      const auto delta = std::sqrt(rij.dot(rij))-(Ri+Rj);  // negative will get ya a force

      // herzian for now
      const auto c1 = 4.0/3.0*mYoungsModulus*std::sqrt(Rij);
      const auto c2 = std::sqrt(4.0*mij*c1/(1+mBeta*mBeta));
      if (delta < 0.0){
        const auto vn = vij.dot(rhatij);
        const auto f = -(c1*delta - c2*vn);
        DvDti += f/mi*rhatij;
        DvDtj -= f/mj*rhatij;
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region


  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
        const auto veli = velocity(nodeListi,i);
        DxDt(nodeListi,i) = veli;
    }
  }

};


}
