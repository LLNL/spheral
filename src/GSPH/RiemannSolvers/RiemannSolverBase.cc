#include "FileIO/FileIO.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Boundary/Boundary.hh"

#include "Hydro/HydroFieldNames.hh"
#include "GSPH/GSPHFieldNames.hh"

#include "GSPH/WaveSpeeds/WaveSpeedBase.hh"
#include "GSPH/Limiters/LimiterBase.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"


#ifdef _OPENMP
#include "omp.h"
#endif

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
RiemannSolverBase<Dimension>::
RiemannSolverBase(LimiterBase<Dimension>& slopeLimiter,
                  WaveSpeedBase<Dimension>& waveSpeed,
                  bool linearReconstruction,
                  int gradType):
  mSlopeLimiter(slopeLimiter),
  mWaveSpeed(waveSpeed),
  mLinearReconstruction(linearReconstruction),
  mGradType(gradType),
  mDpDx(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mDrhoDx(FieldStorageType::CopyFields){

}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
RiemannSolverBase<Dimension>::
~RiemannSolverBase(){}

//========================================================
// initialize derivs -- we'll stowe copies of some spatial 
// derivs. the basic set up will be pressure and velocity
//========================================================
template<typename Dimension>
void
RiemannSolverBase<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           const State<Dimension>& state,
           const StateDerivatives<Dimension>& derivs,
           typename RiemannSolverBase<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename RiemannSolverBase<Dimension>::ConstBoundaryIterator boundaryEnd,
           const typename Dimension::Scalar /*time*/,
           const typename Dimension::Scalar /*dt*/,
           const TableKernel<Dimension>& /*W*/){

  if(mLinearReconstruction){
    dataBase.resizeFluidFieldList(mDpDx,Vector::zero,"riemann solver pressure gradient",true);
    dataBase.resizeFluidFieldList(mDvDx,Tensor::zero,"riemann solver velocity Gradient",true);
    //dataBase.resizeFluidFieldList(mDrhoDx,Vector::zero,"riemann solver density Gradient",true);
    
    //const auto& DrhoDx0 = derivs.fields( GSPHFieldNames::densityGradient, Vector::zero);
    const auto& DpDx0 = derivs.fields( GSPHFieldNames::pressureGradient, Vector::zero);
    const auto& DpDxRaw0 = derivs.fields( GSPHFieldNames::pressureGradient+"RAW", Vector::zero);
    const auto& DvDx0 = derivs.fields( HydroFieldNames::velocityGradient,Tensor::zero);
    const auto& DvDxRaw0 = derivs.fields( HydroFieldNames::velocityGradient+"RAW",Tensor::zero);
    const auto& DvDt0 = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
    const auto& rho0 = state.fields(HydroFieldNames::massDensity, 0.0);

    const auto& connectivityMap = dataBase.connectivityMap();
    const auto& nodeLists = connectivityMap.nodeLists();
    const auto numNodeLists = nodeLists.size();

    // copy from previous time step
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto& nodeList = nodeLists[nodeListi];
      const auto ni = nodeList->numInternalNodes();
      #pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
        const auto DvDxi = DvDx0(nodeListi,i);
        const auto DvDxRawi = DvDxRaw0(nodeListi,i);
        const auto DpDxi = DpDx0(nodeListi,i);
        const auto DpDxRawi = DpDxRaw0(nodeListi,i);
        const auto DvDti = DvDt0(nodeListi,i);
        const auto rhoi = rho0(nodeListi,i);

        // this'll need some cleaning
        switch(mGradType){ 
          case 1 : // default grad based on riemann soln
            mDvDx(nodeListi,i) = DvDxi;
            mDpDx(nodeListi,i) = DpDxi;
            break;
          case 2 : // based on hydro accel for DpDx
            mDvDx(nodeListi,i) = DvDxi;
            mDpDx(nodeListi,i) = -rhoi*DvDti;
            break;
          case 3 : // raw gradients
            mDvDx(nodeListi,i) = DvDxRawi;
            mDpDx(nodeListi,i) = DpDxRawi;
            break;
          case 4 : // raw gradient for P riemann gradient for v
            mDvDx(nodeListi,i) = DvDxi;
            mDpDx(nodeListi,i) = DpDxRawi;
            break;
          case 5 : // raw gradients
            mDvDx(nodeListi,i) = DvDxi;
            mDpDx(nodeListi,i) = Vector::zero;
            break;
          default : 
            mDvDx(nodeListi,i) = Tensor::zero;
            mDpDx(nodeListi,i) = Vector::zero;
            
        }
      } 
    }

    for (auto boundItr = boundaryBegin;
              boundItr != boundaryEnd;
            ++boundItr) {
      (*boundItr)->applyFieldListGhostBoundary(mDpDx);
      (*boundItr)->applyFieldListGhostBoundary(mDvDx);
      //(*boundItr)->applyFieldListGhostBoundary(mDrhoDx);
    }

    for (auto boundItr = boundaryBegin;
              boundItr != boundaryEnd;
            ++boundItr) (*boundItr)->finalizeGhostBoundary();
  
  } // if LinearReconstruction

} // initialize method


//========================================================
// default to non-op
//========================================================
template<typename Dimension>
void
RiemannSolverBase<Dimension>::
interfaceState(const int /*i*/,
               const int /*j*/,
               const int /*nodelisti*/,
               const int /*nodelistj*/,
               const Vector& /*ri*/,
               const Vector& /*rj*/,
               const Scalar& /*rhoi*/,   
               const Scalar& /*rhoj*/, 
               const Scalar& /*ci*/,   
               const Scalar& /*cj*/, 
               const Scalar& /*sigmai*/,    
               const Scalar& /*sigmaj*/,
               const Vector& /*vi*/,    
               const Vector& /*vj*/,
                     Scalar& /*Pstar*/,
                     Vector& /*vstar*/,
                     Scalar& /*rhostari*/,
                     Scalar& /*rhostarj*/) const{

}


//========================================================
// default to non-op
//========================================================
template<typename Dimension>
void
RiemannSolverBase<Dimension>::
interfaceState(const int /*i*/,
               const int /*j*/,
               const int /*nodelisti*/,
               const int /*nodelistj*/,
               const Vector& /*ri*/,
               const Vector& /*rj*/,
               const Scalar& /*rhoi*/,   
               const Scalar& /*rhoj*/, 
               const Scalar& /*ci*/,   
               const Scalar& /*cj*/,
               const Scalar& /*Pi*/,    
               const Scalar& /*Pj*/,
               const Vector& /*vi*/,    
               const Vector& /*vj*/,
               const SymTensor& /*Si*/,    
               const SymTensor& /*Sj*/,
               const Tensor& /*Di*/,    
               const Tensor& /*Dj*/,
                     Vector& /*Tstar*/,
                     Vector& /*vstar*/) const{

}

} // spheral namespace