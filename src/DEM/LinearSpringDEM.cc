//---------------------------------Spheral++----------------------------------//
// DEM -- damped linear spring contact model based on pkdgrav immplementation
//---------------------------------------------------------------------------
// Schwartz, S.R. and Richards, D.C. "An implementation of the soft-sphere 
// discrete element method in a high-performance parallel gravity tree-code,"
// Granular Matter, (2012) 14:363–380, 10.1007/s10035-012-0346-z.
//
// Zhang et. al. "Rotational Failure of Rubble-pile Bodies: Influences of 
// Shear and Cohesive Strengths," The Astrophysical Journal, (2018) 857:15, 20
// 10.3847/1538-4357/aab5b2.
//
//----------------------------------------------------------------------------//
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Hydro/HydroFieldNames.hh"

#include "DEM/ReplaceAndIncrementPairFieldList.hh"
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
  const auto shapeFactor = 1.00;
  const auto muD = 0.3;
  const auto muS = 0.4;
  const auto muT = 1.0;
  const auto muR = 1.0;

  const auto beta2 = shapeFactor*shapeFactor;
  const auto dampingConstTerms = 4.0*mNormalSpringConstant/(1.0+mBeta*mBeta);
  
  // spring constants
  const auto kn = mNormalSpringConstant;           // normal
  const auto ks = 1.0*mNormalSpringConstant;   // sliding
  const auto kt = 2.0*ks*beta2;                    // torsion
  const auto kr = mNormalSpringConstant*beta2;     // rolling
  const auto invKs = 1.0/ks;
  const auto invKt = 1.0/kt;
  const auto invKr = 1.0/kr;
  

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  
  // Get the state and derivative FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto omega = state.fields(DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  const auto radius = state.fields(DEMFieldNames::particleRadius, 0.0);
  const auto equilibriumOverlap = state.fields(DEMFieldNames::equilibriumOverlap, std::vector<Scalar>());
  const auto shearDisplacement = state.fields(DEMFieldNames::shearDisplacement, std::vector<Vector>());
  const auto rollingDisplacement = state.fields(DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  const auto torsionalDisplacement = state.fields(DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  const auto neighborIds = state.fields(DEMFieldNames::neighborIndices, std::vector<int>());
  
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(radius.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);
  CHECK(equilibriumOverlap.size() == numNodeLists);
  CHECK(shearDisplacement.size() == numNodeLists);
  CHECK(rollingDisplacement.size() == numNodeLists);
  CHECK(torsionalDisplacement.size() == numNodeLists);
  CHECK(neighborIds.size() == numNodeLists);

  auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto DomegaDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  auto DDtShearDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() +  DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto newShearDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() +  DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto DDtRollingDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() +  DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto newRollingDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() +  DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto DDtTorsionalDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::incrementPrefix() +  DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  auto newTorsionalDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::replacePrefix() +  DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());

  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DomegaDt.size() == numNodeLists);
  CHECK(DDtShearDisplacement.size() == numNodeLists);
  CHECK(newShearDisplacement.size() == numNodeLists);
  CHECK(DDtRollingDisplacement.size() == numNodeLists);
  CHECK(newRollingDisplacement.size() == numNodeLists);
  CHECK(DDtTorsionalDisplacement.size() == numNodeLists);
  CHECK(newTorsionalDisplacement.size() == numNodeLists);

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
      const auto deltaSlidij = shearDisplacement(storeNodeList,storeNode)[storeContact];
      const auto deltaRollij = rollingDisplacement(storeNodeList,storeNode)[storeContact];
      const auto deltaTorsij = torsionalDisplacement(storeNodeList,storeNode)[storeContact];
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

        // reduced radii
        const auto li = (Ri*Ri-Rj*Rj + rijMag*rijMag)/(2.0*rijMag);
        const auto lj = rijMag-li;

        // effective quantities (mass, reduced radius, contact radius)
        const auto mij = (mi*mj)/(mi+mj);
        const auto lij = 2.0*(li*lj)/(li+lj);
        const auto aij = std::sqrt(4.0*rijMag*Ri - (rijMag*rijMag - Rj*Rj + Ri*Ri))/(2.0*rijMag);

        // damping constants
        const auto Cn = std::sqrt(mij*dampingConstTerms);
        const auto Cs = Cn;
        const auto Ct = 2.0 * Cs * beta2 * lij * lij;
        const auto Cr = Cn * beta2 * lij * lij;

        // our velocities
        const auto vroti = -DEMDimension<Dimension>::cross(omegai,rhatij);
        const auto vrotj =  DEMDimension<Dimension>::cross(omegaj,rhatij);

        const auto vij = vi-vj + li*vroti - lj*vrotj;

        const auto vn = vij.dot(rhatij);
        const auto vt = vij - vn*rhatij;

        // normal force
        const auto fn = (kn*delta - Cn*vn)*rhatij;
        const auto fnMag = fn.magnitude();

        // sliding
        //------------------------------------------------------------
        // project onto new tangential plane -- maintain magnitude
        auto newDeltaSlidij = (deltaSlidij - rhatij.dot(deltaSlidij)*rhatij).unitVector()*deltaSlidij.magnitude();
        
        // spring dashpot
        const auto ft0spring = - ks*newDeltaSlidij;
        const auto ft0damp = - Cs*vt;
        auto ft = ft0spring + ft0damp;

        // static friction limit
        if  (ft.magnitude() > muS*fnMag){

          const auto ftDynamic = muD*fnMag;
          ft = ftDynamic*ft.unitVector();

          newDeltaSlidij = ( ft0damp.magnitude() > ftDynamic ? 
                                Vector::zero : 
                               -invKs*(ft-ft0damp) );

        }

        // Moment - rolling velocity
        //const typename DEMDimension<Dimension>::AngularVector Mroll = -muR*fnMag*this->rollingMoment(rhatij,vroti,vrotj);
        //DomegaDti -= Mroll*lij;
        //DomegaDtj += Mroll*lij;

        // Moment - torsion 
        auto newDeltaTorsij = deltaTorsij;
        const auto omeganij = DEMDimension<Dimension>::dot(omegai-omegaj,rhatij);
        
        // spring dashpot
        const auto Mt0spring = - kt*newDeltaTorsij;
        const auto Mt0damp = - Ct*omeganij;
        auto MtorsionMag = (Mt0spring + Mt0damp)*lij*lij;
        const auto MtStatic = muT*fnMag;

        // limit to static
        if  (MtorsionMag > MtStatic){
          MtorsionMag =  MtStatic;
          newDeltaTorsij = (Mt0damp > MtStatic ? 0.0 :  -invKt*(MtorsionMag-Mt0damp));
        }

        // Rectilinear Acceleration
        const auto fij = fn + ft;
        DvDti += fij;
        DvDtj -= fij;

        // angular acceleration
        const auto Msliding = -DEMDimension<Dimension>::cross(rhatij,fij);
        const auto Mtorsion = MtorsionMag * this->torsionMoment(omegai,omegaj,rhatij);
        DomegaDti += Msliding*li - (Mtorsion)*lij;
        DomegaDtj += Msliding*lj + (Mtorsion)*lij;

        // for spring updates
        newShearDisplacement(storeNodeList,storeNode)[storeContact] = newDeltaSlidij;
        DDtShearDisplacement(storeNodeList,storeNode)[storeContact] = vt;
        newTorsionalDisplacement(storeNodeList,storeNode)[storeContact] = newDeltaTorsij;
        DDtTorsionalDisplacement(storeNodeList,storeNode)[storeContact] = omeganij;
    
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
