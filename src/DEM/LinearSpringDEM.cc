//---------------------------------Spheral++----------------------------------//
// Physics -- dampled linear spring contact model Spheral++
//----------------------------------------------------------------------------//
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DEM/DEMFieldNames.hh"
#include "DEM/DEMDimension.hh"
#include "DEM/LinearSpringDEM.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <cmath>
#include <limits>
namespace Spheral {

//------------------------------------------------------------------------------
// couple explicit things for dimension templating
//------------------------------------------------------------------------------

// moment of interia
// inline
// Dim<1>::Scalar
// momentOfInteria(const Dim<1>::Scalar m, const Dim<1>::Scalar R) {
//   return 1.0;
// }

// inline
// Dim<2>::Scalar
// momentOfInteria(const Dim<2>::Scalar m, const Dim<2>::Scalar R) {
//   return 0.5*m*R*R;
// }

inline
Dim<3>::Scalar
momentOfInteria(const Dim<3>::Scalar m, const Dim<3>::Scalar R) {
  return 0.4*m*R*R;
}

// vel due to rotation
// inline
// Dim<1>::Vector
// crossproduct(const Dim<1>::Vector v1, const Dim<1>::Vector v2) {
//   return Dimension::Vector::zero;
// }

// inline
// Dim<2>::Vector
// crossproduct(const Dim<2>::Vector v1, const Dim<2>::Vector v2) {
//   return Vector();
// }

// inline
// Dim<3>::Vector
// crossproduct(const Dim<3>::Vector v1, const Dim<3>::Vector v2) {
//   return rotation.cross(direction);
// }

//------------------------------------------------------------------------------
// Default constructor
//------------------------------------------------------------------------------
template<typename Dimension>
LinearSpringDEM<Dimension>::
LinearSpringDEM(const DataBase<Dimension>& dataBase,
                   const Scalar normalSpringConstant,
                   const Scalar restitutionCoefficient,
                   const Scalar stepsPerCollision,
                   const TableKernel<Dimension>& W,
                   const Vector& xmin,
                   const Vector& xmax):
                   DEMBase<Dimension>(dataBase,W,stepsPerCollision,xmin,xmax),
                   mNormalSpringConstant(normalSpringConstant),
                   mRestitutionCoefficient(restitutionCoefficient){
     
      const auto pi = 3.14159265358979323846;
      const auto mass = dataBase.DEMMass();
      const auto minMass = mass.min();

      mBeta = pi/std::log(restitutionCoefficient);
      mTimeStep = pi*std::sqrt(0.5*minMass/normalSpringConstant * (1.0 + 1.0/(mBeta*mBeta)));

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
LinearSpringDEM<Dimension>::
~LinearSpringDEM() {}


//------------------------------------------------------------------------------
// time step -- constant for this model
//------------------------------------------------------------------------------
template<typename Dimension>
typename LinearSpringDEM<Dimension>::TimeStepType
LinearSpringDEM<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/,
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const typename Dimension::Scalar /*currentTime*/) const{
  return make_pair(mTimeStep,("DEM vote for time step"));
};

//------------------------------------------------------------------------------
// get our acceleration and other things
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const{

  // A few useful constants we'll use in the following loop.
  //const double tiny = std::numeric_limits<double>::epsilon();
  const auto dampingConstTerms = 4.0*mNormalSpringConstant/(1.0+mBeta*mBeta);

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
  const auto omega = state.fields(DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  const auto radius = state.fields(DEMFieldNames::particleRadius, 0.0);

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(radius.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);

  //auto  T    = derivatives.getany("minContactTime",0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DomegaDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);

  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DomegaDt.size() == numNodeLists);


#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    //auto DomegaDt_thread = DomegaDt.threadCopy(threadStack);

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
      const auto& omegai = omega(nodeListi, i);
      const auto& Ri = radius(nodeListi, i);
      
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DomegaDti = DomegaDt(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& omegaj = omega(nodeListj, j);
      const auto& Rj = radius(nodeListj, j);

      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DomegaDtj = DomegaDt(nodeListj, j);

      CHECK(mi > 0.0);
      CHECK(mj > 0.0);
      CHECK(Ri > 0.0);
      CHECK(Rj > 0.0);

      // are we overlapping ? 
      const auto rij = ri-rj;
      const auto rijMag = rij.magnitude();
      const auto delta = (Ri+Rj) - rijMag; 
      
      // if so do the things
      if (delta > 0.0){
      
        // line of action for the contact
        const auto rhatij = rij.unitVector();

        //total relative velocity
        //const auto li = (Ri*Ri-Rj*Rj + rijMag*rijMag)/(2.0*rijMag);
        //const auto lj = rijMag-li;

        const auto vij = vi-vj;//+ lj*() - li*(rhatij.cross(omegai));
        const auto vn = vij.dot(rhatij);
        const auto vt = vij - vn*rhatij;

        // effective quantities
        const auto mij = (mi*mj)/(mi+mj);
        
        // moments of interia
        const auto Ii = momentOfInteria(mi,Ri);
        const auto Ij = momentOfInteria(mj,Rj);

        // normal force w/ Herzian spring constant
        const auto normalDampingConstant = std::sqrt(mij*dampingConstTerms);

        // normal force
        const auto f = mNormalSpringConstant*delta - normalDampingConstant*vn;
        //const auto M = 

        DvDti += f/mi*rhatij;
        DvDtj -= f/mj*rhatij;

        //DomegaDti += M/I;
      }  
    } // loop over pairs
    threadReduceFieldLists<Dimension>(threadStack);
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
