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
#include "DEM/DEMFieldNames.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <cmath>
#include <limits>
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


//------------------------------------------------------------------------------
// time step 
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
DampedLinearSpring<Dimension>::
timeStep(const DataBase<Dimension>& dataBase,
         const State<Dimension>& state,
         const StateDerivatives<Dimension>& /*derivs*/,
         const typename Dimension::Scalar /*currentTime*/) const{
  
  const auto& mass = state.fields(HydroFieldNames::mass, 0.0); 
  const auto& radius = state.fields(DEMFieldNames::particleRadius,0.0);

  const auto numNodeLists = dataBase.numNodeLists();

  const auto pi = 3.14159265358979323846;
  auto minContactTime = std::numeric_limits<double>::max();

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
  const auto& nodeList = mass[nodeListi]->nodeList();
  const auto ni = nodeList.numInternalNodes();

  #pragma omp parallel
  {
    auto minContactTime_thread = std::numeric_limits<double>::max();
    #pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
          const auto mi = mass(nodeListi,i);
          const auto Ri = radius(nodeListi,i);
          const auto k = 4.0/3.0*mYoungsModulus*std::sqrt(Ri);
          minContactTime_thread = min(minContactTime_thread,mi/k);
      }

    #pragma omp critical
    if (minContactTime_thread < minContactTime) minContactTime = minContactTime_thread;

    }
  }

  minContactTime = pi*std::sqrt(minContactTime);
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
  const double tiny = std::numeric_limits<double>::epsilon();
  //auto minContactTime = std::numeric_limits<double>::max();

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
  //const auto omega = state.fields(DEMFieldNames::angularVelocity, Vector::zero);
  const auto radius = state.fields(DEMFieldNames::particleRadius, .0);

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(radius.size() == numNodeLists);
  //CHECK(omega.size() == numNodeLists);

  //auto  T    = derivatives.getany("minContactTime",0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  //auto  DomegaDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + DEMFieldNames::angularVelocity, Vector::zero);

  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  //CHECK(DomegaDt.size() == numNodeLists);


#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    //auto minContactTimekk = std::numeric_limits<double>::max();

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
      CHECK(Ri > 0.0);
      CHECK(mj > 0.0);
      CHECK(Hdetj > 0.0);
      CHECK(Rj > 0.0);
      
      const auto mij = (mi*mj)/(mi+mj);
      const auto Rij = (Ri*Rj)/(Ri+Rj);

      const auto vij = vi-vj;
      const auto rij = ri-rj;
      const auto rhatij = rij.unitVector();

      const auto delta = (Ri+Rj) - std::sqrt(rij.dot(rij));  // negative will get ya a force
      
      const auto c1 = 4.0/3.0*mYoungsModulus*std::sqrt(Rij);
      const auto c2 = std::sqrt(4.0*mij*c1/(1.0+mBeta*mBeta));
      
      if (delta > 0.0){
         const auto vn = vij.dot(rhatij);
         const auto f = c1*delta - c2*vn;
         DvDti += f/mi*rhatij;
         DvDtj -= f/mj*rhatij;
      }

      //minContactTimekk = min(c1/mij,minContactTimekk);
      
    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
    //#pragma omp critical
    //  if (minContactTimekk < minContactTime) minContactTime = minContactTimekk;
  }   // OpenMP parallel region

  
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
        const auto veli = velocity(nodeListi,i);
        DxDt(nodeListi,i) = veli;
    }
  }

};


}
