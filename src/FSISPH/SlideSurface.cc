//---------------------------------Spheral++----------------------------------//
// Slide Surface -- 
//----------------------------------------------------------------------------//
#include "FSISPH/SlideSurface.hh"
#include "FSISPH/computeSurfaceNormals.hh"
#include "FSISPH/computeSurfaceSmoothness.hh"
#include "FSISPH/FSIFieldNames.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"

#include "FileIO/FileIO.hh"

#include "Hydro/HydroFieldNames.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"

#include "Boundary/Boundary.hh"

#include "Neighbor/ConnectivityMap.hh"

#include <vector>
#include <iostream>
namespace Spheral {

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------
template<typename Dimension>
SlideSurface<Dimension>::
SlideSurface(const TableKernel<Dimension>& W,
             const std::vector<int> contactTypes):
  Physics<Dimension>(),
  mKernel(W),
  mIsActive(false),
  mNumNodeLists(0.0),
  mIsSlideSurface(),
  mSurfaceNormals(FieldStorageType::CopyFields),
  mSurfaceFraction(FieldStorageType::CopyFields),
  mSurfaceSmoothness(FieldStorageType::CopyFields){

    for(std::vector<int>::const_iterator it = contactTypes.begin();
        it != contactTypes.end();
        ++it){
      if (*it == 1){
        mIsActive=true;
        mIsSlideSurface.push_back(true);
      }else{
        mIsSlideSurface.push_back(false);
      }    
    }

}

template<typename Dimension>
SlideSurface<Dimension>::
SlideSurface(const TableKernel<Dimension>& W,
             const std::vector<bool> contactTypes):
  Physics<Dimension>(),
  mKernel(W),
  mIsActive(false),
  mNumNodeLists(0.0),
  mIsSlideSurface(contactTypes),
  mSurfaceNormals(FieldStorageType::CopyFields),
  mSurfaceFraction(FieldStorageType::CopyFields),
  mSurfaceSmoothness(FieldStorageType::CopyFields){

    for(std::vector<bool>::const_iterator it = contactTypes.begin();
        it != contactTypes.end();
        ++it){
      if (*it) mIsActive=true;  
    }

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SlideSurface<Dimension>::
~SlideSurface() {
}

//------------------------------------------------------------------------------
// more intelligable access
//------------------------------------------------------------------------------
template<typename Dimension>
bool 
SlideSurface<Dimension>::
isSlideSurface(const int nodeListi, 
               const int nodeListj) const {
    const auto oneDimIndex = mNumNodeLists * nodeListi + nodeListj;
    return mIsSlideSurface[oneDimIndex];
};

//------------------------------------------------------------------------------
// return correction factor [0,1] for the artificial viscosity pressure
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar 
SlideSurface<Dimension>::
slideCorrection(const int nodeListi,
                const int i, 
                const int nodeListj,
                const int j,
                const typename Dimension::Vector vi,
                const typename Dimension::Vector vj) const {

    typename Dimension::Scalar slideCorr = 1.0;
    if (this->isSlideSurface(nodeListi,nodeListj)){
      
      const auto normi = mSurfaceNormals(nodeListi,i);
      const auto normj = mSurfaceNormals(nodeListj,j);

      // surface smoothness to turn of our slide 
      // 1.00-0.95 full on
      // 0.95-0.70 linear ramp down
      // 0.70-0.00 full off
      const auto maxSS = min(mSurfaceSmoothness(nodeListi,i),mSurfaceSmoothness(nodeListj,j));
      const auto ssij  = 4.0*min(max(0.95-maxSS,0.0),0.25);

      const auto vijhat = (vi-vj).unitVector();
      const auto fi = abs(normi.dot(vijhat));//max(0.0, normi.dot(vijhat));
      const auto fj = abs(normj.dot(vijhat));//max(0.0,-normj.dot(vijhat));
      slideCorr = ssij + (1.0-ssij)*fi*fj; 
    
    }
    // if (this->isSlideSurface(nodeListi,nodeListj)){
    //   const auto normi = mSurfaceNormals(nodeListi,i);
    //   const auto normj = mSurfaceNormals(nodeListj,j);
    //   const auto vijhat = (vi-vj).unitVector();
    //   const auto fi =  abs(normi.dot(vijhat));
    //   const auto fj =  abs(normj.dot(vijhat));
    //   slideCorr = fi*fj; 
    // }
    return slideCorr;      

}

//------------------------------------------------------------------------------
// called once on problem start
//------------------------------------------------------------------------------
template<typename Dimension>
void
SlideSurface<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase){
  mSurfaceNormals = dataBase.newFluidFieldList(Vector::zero,  FSIFieldNames::interfaceNormals);
  mSurfaceFraction = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceFraction);
  mSurfaceSmoothness = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceSmoothness);
  mNumNodeLists = dataBase.numNodeLists();
}

//------------------------------------------------------------------------------
// normals calculated prior to eval derivs step
//------------------------------------------------------------------------------
template<typename Dimension>
void
SlideSurface<Dimension>::
initialize(const typename Dimension::Scalar /*time*/,
           const typename Dimension::Scalar /*dt*/,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& /*derivs*/) {

const auto& connectivityMap = dataBase.connectivityMap();
  const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto& massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
        auto  normals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
        auto  surfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
        auto  surfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

  surfaceFraction.Zero();
  surfaceSmoothness.Zero();
  normals.Zero();
  
  computeSurfaceNormals(connectivityMap,
                        mKernel,
                        position,
                        mass,
                        massDensity,
                        H,
                        normals);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr)(*boundaryItr)->applyFieldListGhostBoundary(normals);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  
  computeSurfaceSmoothness(connectivityMap,
                             mKernel,
                             position,
                             mass,
                             massDensity,
                             H,
                             normals,
                             surfaceFraction,
                             surfaceSmoothness);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr){
    (*boundaryItr)->applyFieldListGhostBoundary(surfaceSmoothness);
    (*boundaryItr)->applyFieldListGhostBoundary(surfaceFraction);
  }


  
};


//------------------------------------------------------------------------------
// register the surface normals w/ the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
SlideSurface<Dimension>:: 
registerState(DataBase<Dimension>& dataBase,
                   State<Dimension>& state){
  dataBase.resizeFluidFieldList(mSurfaceNormals, Vector::zero, FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mSurfaceFraction, 0.0, FSIFieldNames::interfaceFraction,false); 
  dataBase.resizeFluidFieldList(mSurfaceSmoothness, 0.0, FSIFieldNames::interfaceSmoothness,false);
  state.enroll(mSurfaceNormals); 
  state.enroll(mSurfaceFraction);
  state.enroll(mSurfaceSmoothness);               
};


//------------------------------------------------------------------------------
// non-op methods from Physics
//------------------------------------------------------------------------------
template<typename Dimension> 
typename SlideSurface<Dimension>::TimeStepType  
SlideSurface<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const typename Dimension::Scalar /*currentTime*/) const{

   return make_pair(std::numeric_limits<double>::max(), this->label());
};



template<typename Dimension>
void
SlideSurface<Dimension>:: 
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                          StateDerivatives<Dimension>& /*derivatives*/) const{};

template<typename Dimension>
void
SlideSurface<Dimension>:: 
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& /*derivs*/){};


template<typename Dimension>
void
SlideSurface<Dimension>::
dumpState(FileIO& /*file*/, 
          const std::string& /*pathName*/) const{};

template<typename Dimension>
void
SlideSurface<Dimension>::
restoreState(const FileIO& /*file*/,
             const std::string& /*pathName*/){};

}