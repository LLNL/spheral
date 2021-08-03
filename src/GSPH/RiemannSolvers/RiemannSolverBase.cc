#include "FileIO/FileIO.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"

#include "Hydro/HydroFieldNames.hh"
#include "GSPH/GSPHFieldNames.hh"

#include "GSPH/WaveSpeeds/WaveSpeedBase.hh"
#include "GSPH/Limiters/SlopeLimiterBase.hh"
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
RiemannSolverBase(SlopeLimiterBase<Dimension>& slopeLimiter,
                  WaveSpeedBase<Dimension>& waveSpeed):
  mSlopeLimiter(slopeLimiter),
  mWaveSpeed(waveSpeed),
  mDpDx(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields){

}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
RiemannSolverBase<Dimension>::
~RiemannSolverBase(){}

//========================================================
// initialize derivs -- we'll stowe copies of some spatial 
// derivs
//========================================================
template<typename Dimension>
void
RiemannSolverBase<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
                const State<Dimension>& state,
                const StateDerivatives<Dimension>& derivs,
                const typename Dimension::Scalar time,
                const typename Dimension::Scalar dt,
                const TableKernel<Dimension>& W){

  dataBase.resizeFluidFieldList(mDpDx,Vector::zero,"pressureGradient",true);
  dataBase.resizeFluidFieldList(mDvDx,Tensor::zero,"velocityGradient",true);

  const auto DpDx = derivs.fields( GSPHFieldNames::pressureGradient, Vector::zero);
  const auto DvDx = derivs.fields( HydroFieldNames::velocityGradient,Tensor::zero);

  
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // copy from previous time step
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = nodeLists[nodeListi];
    const auto ni = nodeList->numNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto DvDxi = DvDx(nodeListi,i);
      const auto DpDxi = DpDx(nodeListi,i);
      mDvDx(nodeListi,i) = DvDxi;
      mDpDx(nodeListi,i) = DpDxi;
    } 
  }

} // initialize method 

} // spheral namespace