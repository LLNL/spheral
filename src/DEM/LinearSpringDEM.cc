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
#include "FileIO/FileIO.hh"

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

#include <iostream>
#include <stdexcept>
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
                const Scalar normalRestitutionCoefficient,
                const Scalar tangentialSpringConstant,
                const Scalar tangentialRestitutionCoefficient,
                const Scalar dynamicFrictionCoefficient,
                const Scalar staticFrictionCoefficient,
                const Scalar rollingFrictionCoefficient,
                const Scalar torsionalFrictionCoefficient,
                const Scalar cohesiveTensileStrength,
                const Scalar shapeFactor,
                const Scalar stepsPerCollision,
                const Vector& xmin,
                const Vector& xmax):
  DEMBase<Dimension>(dataBase,stepsPerCollision,xmin,xmax),
  mNormalSpringConstant(normalSpringConstant),
  mNormalRestitutionCoefficient(normalRestitutionCoefficient),
  mTangentialSpringConstant(tangentialSpringConstant),
  mTangentialRestitutionCoefficient(tangentialRestitutionCoefficient),
  mDynamicFrictionCoefficient(dynamicFrictionCoefficient),
  mStaticFrictionCoefficient(staticFrictionCoefficient),
  mRollingFrictionCoefficient(rollingFrictionCoefficient),
  mTorsionalFrictionCoefficient(torsionalFrictionCoefficient),
  mCohesiveTensileStrength(cohesiveTensileStrength),
  mShapeFactor(shapeFactor){

  this->setTimeStep(dataBase);

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
LinearSpringDEM<Dimension>::
~LinearSpringDEM() {}

//------------------------------------------------------------------------------
// set our timestep (for now its a constant single value)
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
setTimeStep(const DataBase<Dimension>& dataBase){
    
    const auto pi = 3.14159265358979323846;
    const auto mass = dataBase.DEMMass();
    const auto minMass = mass.min();

    mNormalBeta = pi/std::log(std::max(mNormalRestitutionCoefficient,1.0e-3));
    mTangentialBeta = pi/std::log(std::max(mTangentialRestitutionCoefficient,1.0e-3));
    mTimeStep = pi*std::sqrt(0.5*minMass/mNormalSpringConstant * (1.0 + 1.0/(mNormalBeta*mNormalBeta)));

}


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
}

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
  const double tiny = std::numeric_limits<double>::epsilon();
  const auto shapeFactor = this->shapeFactor();
  const auto shapeFactor2 = shapeFactor*shapeFactor; 

  const auto muD = this->dynamicFrictionCoefficient();
  const auto muS = this->staticFrictionCoefficient();
  const auto muT = this->torsionalFrictionCoefficient();
  const auto muR = this->rollingFrictionCoefficient();

  const auto Cc = this->cohesiveTensileStrength();

  // spring constants
  const auto kn = this->normalSpringConstant();        // normal
  const auto ks = this->tangentialSpringConstant();    // sliding
  const auto kt = 0.50 * ks * shapeFactor2;
  const auto kr = 0.25 * kn * shapeFactor2;

  const auto invKs = 1.0/max(ks,tiny);
  const auto invKt = 1.0/max(kt,tiny);
  const auto invKr = 1.0/max(kr,tiny);
  
  const auto normalDampingTerms = 2.0*kn/(1.0+mNormalBeta*mNormalBeta);
  const auto tangentialDampingTerms = 4.0/5.0*ks/(1.0+mTangentialBeta*mTangentialBeta);
 
  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  
  // Get the state FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto omega = state.fields(DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  const auto radius = state.fields(DEMFieldNames::particleRadius, 0.0);
  const auto uniqueIndices = state.fields(DEMFieldNames::uniqueIndices, (int)0);
  
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(radius.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);
  CHECK(uniqueIndices.size() == numNodeLists);

  // Get the state pairFieldLists.
  const auto equilibriumOverlap = state.fields(DEMFieldNames::equilibriumOverlap, std::vector<Scalar>());
  const auto shearDisplacement = state.fields(DEMFieldNames::shearDisplacement, std::vector<Vector>());
  const auto rollingDisplacement = state.fields(DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  const auto torsionalDisplacement = state.fields(DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  const auto neighborIds = state.fields(DEMFieldNames::neighborIndices, std::vector<int>());

  CHECK(equilibriumOverlap.size() == numNodeLists);
  CHECK(shearDisplacement.size() == numNodeLists);
  CHECK(rollingDisplacement.size() == numNodeLists);
  CHECK(torsionalDisplacement.size() == numNodeLists);
  CHECK(neighborIds.size() == numNodeLists);

  // Get the deriv FieldLists
  auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto DomegaDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DomegaDt.size() == numNodeLists);
  
  // Get the deriv pairFieldLists
  auto DDtShearDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() +  DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto newShearDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() +  DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto DDtRollingDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() +  DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto newRollingDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() +  DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto DDtTorsionalDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::incrementPrefix() +  DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  auto newTorsionalDisplacement = derivatives.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::replacePrefix() +  DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());

  CHECK(DDtShearDisplacement.size() == numNodeLists);
  CHECK(newShearDisplacement.size() == numNodeLists);
  CHECK(DDtRollingDisplacement.size() == numNodeLists);
  CHECK(newRollingDisplacement.size() == numNodeLists);
  CHECK(DDtTorsionalDisplacement.size() == numNodeLists);
  CHECK(newTorsionalDisplacement.size() == numNodeLists);

  // storage locations in our pairwise FieldLists
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
      const int pairNodeList = contacts[kk].pairNodeList;
      const int pairNode = contacts[kk].pairNode;

      // simple test for our store/pair selection
      const bool orientationOne = (storeNodeList==nodeListi and storeNode==i and pairNodeList==nodeListj and pairNode==j);
      const bool orientationTwo = (storeNodeList==nodeListj and storeNode==j and pairNodeList==nodeListi and pairNode==i);
      if(not(orientationOne or orientationTwo)) throw std::invalid_argument("AddPositiveIntegers arguments must be positive");

      
      // storage sign to make pairwise values i-j independent
      const auto uIds = uniqueIndices(storeNodeList,storeNode);
      const auto uIdp = uniqueIndices(pairNodeList,pairNode);
      const int storageSign = (uIds <= uIdp ? 1 : -1);

      // stored pair-wise values
      const auto overlapij = equilibriumOverlap(storeNodeList,storeNode)[storeContact];
      const auto deltaSlidij = storageSign*shearDisplacement(storeNodeList,storeNode)[storeContact];
      const auto deltaRollij =             rollingDisplacement(storeNodeList,storeNode)[storeContact];
      const auto deltaTorsij = storageSign*torsionalDisplacement(storeNodeList,storeNode)[storeContact];

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

        // effective quantities (mass, reduced radius) respectively, mij->mi for mi=mj
        const auto mij = 2.0*(mi*mj)/(mi+mj);
        const auto lij = 2.0*(li*lj)/(li+lj);

        // damping constants -- ct and cr derived quantities ala Zhang 2017
        const auto Cn = std::sqrt(mij*normalDampingTerms);
        const auto Cs = std::sqrt(mij*tangentialDampingTerms);
        const auto Ct = 0.50 * Cs * shapeFactor2;
        const auto Cr = 0.25 * Cn * shapeFactor2;

        // our velocities
        const Vector vroti = -DEMDimension<Dimension>::cross(omegai,rhatij);
        const Vector vrotj =  DEMDimension<Dimension>::cross(omegaj,rhatij);

        const Vector vij = vi-vj + li*vroti - lj*vrotj;

        const Scalar vn = vij.dot(rhatij);                  // normal velocity
        const Vector vs = vij - vn*rhatij;                  // sliding velocity
        const Vector vr = -li*vroti - lj*vrotj;              // rolling velocity
        const Scalar vt = -lij*DEMDimension<Dimension>::dot(omegai-omegaj,rhatij); // torsion velocity

        // normal forces 
        //------------------------------------------------------------
        const Vector fn = (kn*delta - Cn*vn)*rhatij;        // normal spring
        const Vector fc = Cc*shapeFactor2*lij*lij*rhatij;   // normal cohesion
        const Scalar fnMag = fn.magnitude();                // magnitude of normal spring force

        // sliding
        //------------------------------------------------------------
        // project onto new tangential plane -- maintain magnitude
        Vector newDeltaSlidij = (deltaSlidij - rhatij.dot(deltaSlidij)*rhatij).unitVector()*deltaSlidij.magnitude();
        
        // spring dashpot
        const Vector ft0spring = - ks*newDeltaSlidij;
        const Vector ft0damp = - Cs*vs;
              Vector ft = ft0spring + ft0damp;

        // static friction limit
        if  (ft.magnitude() > muS*fnMag){

          const Scalar ftDynamic = muD*fnMag;
          ft = ftDynamic*ft.unitVector();

          newDeltaSlidij = ( ft0damp.magnitude() > ftDynamic ? 
                                Vector::zero : 
                               -(ft-ft0damp)*invKs );

        }

        // torsion
        //------------------------------------------------------------
        // since we use a scalar no need to modify here
        auto newDeltaTorsij = deltaTorsij;
        
        // spring dashpot
        const Scalar Mt0spring = - kt*newDeltaTorsij;
        const Scalar Mt0damp = - Ct*vt;
              Scalar MtorsionMag = (Mt0spring + Mt0damp); 
        const Scalar MtStatic = muT*shapeFactor*muS*fnMag;

        // limit to static
        if  (std::abs(MtorsionMag) > MtStatic){
          MtorsionMag =  (MtorsionMag > 0.0 ? 1.0 : -1.0)*MtStatic;
          newDeltaTorsij = (std::abs(Mt0damp) > MtStatic ? 0.0 :  -(MtorsionMag-Mt0damp)*invKt);
        }

        // rolling
        //------------------------------------------------------------
        // project onto new tangential plane -- maintain magnitude
        Vector newDeltaRollij = (deltaRollij - rhatij.dot(deltaRollij)*rhatij).unitVector()*deltaRollij.magnitude();
    
        // spring dashpot
        const Vector Mr0spring = - kr*newDeltaRollij;
        const Vector Mr0damp = - Cr*vr;
              Vector effectiveRollingForce = (Mr0spring + Mr0damp); 
        const Scalar MrStatic = muR*shapeFactor*fnMag;

        // limit to static
        if  (effectiveRollingForce.magnitude() > MrStatic){
         effectiveRollingForce =  MrStatic*effectiveRollingForce.unitVector();
         newDeltaRollij = (Mr0damp.magnitude() > MrStatic ?
                           Vector::zero :  
                           -(effectiveRollingForce-Mr0damp)*invKr);
        }

        // accelerations
        //------------------------------------------------------------
        // Rectilinear Acceleration 
        const Vector fij = fn - fc + ft;
        DvDti += fij;
        DvDtj -= fij;

        // angular acceleration
        const auto Msliding = -DEMDimension<Dimension>::cross(rhatij,fij);
        const auto Mrolling = -DEMDimension<Dimension>::cross(rhatij,effectiveRollingForce);
        const auto Mtorsion = MtorsionMag * this->torsionMoment(rhatij,omegai,omegaj); // rename torsionDirection
        DomegaDti += Msliding*li - (Mtorsion + Mrolling) * lij;
        DomegaDtj += Msliding*lj + (Mtorsion + Mrolling) * lij;

        // for spring updates
        newShearDisplacement(storeNodeList,storeNode)[storeContact] = storageSign*newDeltaSlidij;
        DDtShearDisplacement(storeNodeList,storeNode)[storeContact] = storageSign*vs;
        newTorsionalDisplacement(storeNodeList,storeNode)[storeContact] = storageSign*newDeltaTorsij;
        DDtTorsionalDisplacement(storeNodeList,storeNode)[storeContact] = storageSign*vt;
        newRollingDisplacement(storeNodeList,storeNode)[storeContact] = newDeltaRollij;
        DDtRollingDisplacement(storeNodeList,storeNode)[storeContact] = vr;
    
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

}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  DEMBase<Dimension>::dumpState(file, pathName);
  file.write(mTimeStep, pathName + "/timeStep");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSpringDEM<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  DEMBase<Dimension>::restoreState(file, pathName);
  file.read(mTimeStep, pathName + "/timeStep");
}

}

