//---------------------------------Spheral++----------------------------------//
// Physics -- dampled linear spring contact model Spheral++
//----------------------------------------------------------------------------//
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Hydro/HydroFieldNames.hh"

#include "DEM/ReplaceAndIncrementFieldList.hh"
#include "DEM/DEMFieldNames.hh"
#include "DEM/DEMDimension.hh"
#include "DEM/LinearSpringDEM.hh"
#include "DEM/ContactStorageLocation.hh"

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
LinearSpringDEM<Dimension>::
LinearSpringDEM(const DataBase<Dimension>& dataBase,
                const Scalar normalSpringConstant,
                const Scalar restitutionCoefficient,
                const Scalar stepsPerCollision,
                const Vector& xmin,
                const Vector& xmax):
  DEMBase<Dimension>(dataBase,stepsPerCollision,xmin,xmax),
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
  return make_pair(mTimeStep/this->stepsPerCollision(),("Linear Spring DEM vote for time step"));
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

  this->resizeDerivativePairFieldLists(derivatives);

  // A few useful constants we'll use in the following loop.
  //const double tiny = std::numeric_limits<double>::epsilon();
  const auto dampingConstTerms = 4.0*mNormalSpringConstant/(1.0+mBeta*mBeta);
  const auto tangentialSpringConstant = mNormalSpringConstant;
  const auto invTangentialSpringConstant = 1.0/tangentialSpringConstant;
  const auto muD = 0.3;
  const auto muS = 0.4;
  const auto muT = 0.25;
  const auto muR = 0.1;

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
  const auto equilibriumOverlap = state.fields(DEMFieldNames::equilibriumOverlap, std::vector<Scalar>());
  const auto shearDisplacement = state.fields(DEMFieldNames::shearDisplacement, std::vector<Vector>());
  const auto neighborIds = state.fields(DEMFieldNames::neighborIndices, std::vector<int>());
  
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(radius.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);

  auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto DomegaDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  auto DDtShearDisplacement = derivatives.fields(ReplaceAndIncrementFieldList<Dimension, std::vector<Vector>>::incrementPrefix() +  DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto newShearDisplacement = derivatives.fields(ReplaceAndIncrementFieldList<Dimension, std::vector<Vector>>::replacePrefix() +  DEMFieldNames::shearDisplacement, std::vector<Vector>());

  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DomegaDt.size() == numNodeLists);
  CHECK(DDtShearDisplacement.size() == numNodeLists);
  CHECK(newShearDisplacement.size() == numNodeLists);

  const auto& contacts = this->contactStorageIndices();

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DomegaDt_thread = DomegaDt.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {

      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;
      
      const int storeNodeList = contacts[kk].storeNodeList;
      const int storeNode = contacts[kk].storeNode;
      const int storeContact = contacts[kk].storeContact;

      // stored pair-wise values
      const auto overlapij = equilibriumOverlap(storeNodeList,storeNode)[storeContact];
      const auto sij = shearDisplacement(storeNodeList,storeNode)[storeContact];
      //const int numContacts = neighborIds(storeNodeList,storeNode).size();
      
      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& omegai = omega(nodeListi, i);
      const auto& Ri = radius(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& omegaj = omega(nodeListj, j);
      const auto& Rj = radius(nodeListj, j);

      // Get the derivs from node i
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DomegaDti = DomegaDt_thread(nodeListi, i);

      // Get the derivs from node j
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DomegaDtj = DomegaDt_thread(nodeListj, j);

      CHECK(mi > 0.0);
      CHECK(mj > 0.0);
      CHECK(Ri > 0.0);
      CHECK(Rj > 0.0);

      // are we overlapping ? 
      const auto rij = ri-rj;
      const auto rijMag = rij.magnitude();
      const auto delta0 = (Ri+Rj) - rijMag;  // raw delta
      
      // if so do the things
      if (delta0 > 0.0){

        // effective delta
        const auto delta = delta0 - overlapij;

        // line of action for the contact
        const auto rhatij = rij.unitVector();

        // effective radii
        const auto li = (Ri*Ri-Rj*Rj + rijMag*rijMag)/(2.0*rijMag);
        const auto lj = rijMag-li;

        // effective quantities
        const auto mij = (mi*mj)/(mi+mj);
        const auto lij = (li*lj)/(li+lj);

        // this is super ugly right now and should change down the line...
        const auto vroti = -DEMDimension<Dimension>::cross(omegai,rhatij);
        const auto vrotj =  DEMDimension<Dimension>::cross(omegaj,rhatij);

        const auto vij = vi-vj + li*vroti - lj*vrotj;

        const auto vn = vij.dot(rhatij);
        const auto vt = vij - vn*rhatij;

        // projection of shear spring displacement
        auto newSij = sij - rhatij.dot(sij)*rhatij;
        const auto shatij = newSij.unitVector();
        newSij = shatij * sij.magnitude();

        // normal force coefficient
        const auto normalDampingConstant = std::sqrt(mij*dampingConstTerms);
        const auto tangentialDampingConstant = normalDampingConstant;

        // normal force
        const auto fn = (mNormalSpringConstant*delta - normalDampingConstant*vn)*rhatij;
        const auto fnMag = fn.magnitude();

        // friction force
        const auto ft0damp = - tangentialDampingConstant*vt;
        const auto ft0spring = - tangentialSpringConstant*newSij;
        const auto ftStatic = muS*fnMag;

        auto ft = ft0spring + ft0damp;

        if  (ft.magnitude() > ftStatic){

          const auto ftDynamic = muD*fnMag;
          ft = ftDynamic*ft.unitVector();

          newSij = ( ft0damp.magnitude() > ftDynamic ? 
                     Vector::zero : 
                     -invTangentialSpringConstant*(ft-ft0damp) );

        }

        // pairwise force and moment
        const auto fij = fn + ft;

        // our accelerations
        DvDti += fij;
        DvDtj -= fij;

        // Moment - tangential forces
        const auto Mfriction = -DEMDimension<Dimension>::cross(rhatij,fij);
        DomegaDti += Mfriction*li;
        DomegaDtj += Mfriction*lj;

        // Moment - rolling velocity
        const typename DEMDimension<Dimension>::AngularVector Mroll = -muR*fnMag*this->rollingMoment(rhatij,vroti,vrotj);
        DomegaDti -= Mroll*lij;
        DomegaDtj += Mroll*lij;

        // Moment - torsion 
        const auto contactRadiusij = std::sqrt(4.0*rijMag*Ri - (rijMag*rijMag - Rj*Rj + Ri*Ri))/(2.0*rijMag);
        const auto Mtorsion = muT*fnMag*this->torsionMoment(rhatij,omegai,omegaj);
        DomegaDti -= (Mtorsion)*contactRadiusij;
        DomegaDtj += (Mtorsion)*contactRadiusij;

        newShearDisplacement(storeNodeList,storeNode)[storeContact]=newSij;
        DDtShearDisplacement(storeNodeList,storeNode)[storeContact]=vt;
    
      }  
    } // loop over pairs
    threadReduceFieldLists<Dimension>(threadStack);
  }   // OpenMP parallel region

  
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
        const auto mi = mass(nodeListi,i);
        const auto Ri = radius(nodeListi,i);
        const auto veli = velocity(nodeListi,i);
        const auto Ii = this->momentOfInertia(mi,Ri);

        DxDt(nodeListi,i) = veli;
        DomegaDt(nodeListi,i) /= Ii;
        DvDt(nodeListi,i) /= mi;
    }
  }

};


}
