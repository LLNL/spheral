//---------------------------------Spheral++----------------------------------//
// Slide Surface -- 
//----------------------------------------------------------------------------//

#include "FSISPH/SlideSurface.hh"
#include "FSISPH/FSIFieldNames.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"

#include "FileIO/FileIO.hh"

#include "Hydro/HydroFieldNames.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"

#include "Boundary/Boundary.hh"

#include "Neighbor/ConnectivityMap.hh"

#include "Kernel/TableKernel.hh"

#include <limits.h>
#include <vector>
#include <iostream>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors
//------------------------------------------------------------------------------
template<typename Dimension>
SlideSurface<Dimension>::
SlideSurface(DataBase<Dimension>& dataBase,
             const std::vector<int> contactTypes, 
             const SurfaceNormalMethod surfaceNormalMethod,
             const bool normalsAreSmoothed,
             const bool gradientsAreCorrected):
  mIsActive(false),
  mNormalsAreSmoothed(normalsAreSmoothed),
  mGradientsAreCorrected(gradientsAreCorrected),
  mNumNodeLists(0.0),
  mIsSlideSurface(){

    mNumNodeLists = dataBase.numNodeLists();

    // for our custom "map" (nodelisti,nodelistj) -> bool isSlide
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

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SlideSurface<Dimension>::
~SlideSurface() {
}


//------------------------------------------------------------------------------
// register the surface normals w/ the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
SlideSurface<Dimension>:: 
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state){
  // dataBase.resizeFluidFieldList(mSurfaceNormals, Vector::zero, FSIFieldNames::interfaceNormals,false);
  // dataBase.resizeFluidFieldList(mSurfaceFraction, 0.0, FSIFieldNames::interfaceFraction,false); 
  // dataBase.resizeFluidFieldList(mSurfaceSmoothness, 0.0, FSIFieldNames::interfaceSmoothness,false);
  // state.enroll(mSurfaceNormals); 
  // state.enroll(mSurfaceFraction);
  // state.enroll(mSurfaceSmoothness);               
};

//------------------------------------------------------------------------------
// register the derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
SlideSurface<Dimension>:: 
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs){
  // dataBase.resizeFluidFieldList(mNewSurfaceNormals, Vector::zero,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals,false);
  // dataBase.resizeFluidFieldList(mNewSurfaceFraction, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() +  FSIFieldNames::interfaceFraction,false); 
  // dataBase.resizeFluidFieldList(mNewSurfaceSmoothness, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness,false);
  // dataBase.resizeFluidFieldList(mSmoothedSurfaceNormals, Vector::zero,  FSIFieldNames::smoothedInterfaceNormals,false);
  // dataBase.resizeFluidFieldList(mSmoothnessNormalization, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::smoothnessNormalization,false);
  
  // derivs.enroll(mNewSurfaceNormals); 
  // derivs.enroll(mNewSurfaceFraction);
  // derivs.enroll(mNewSurfaceSmoothness);  
  // derivs.enroll(mSmoothedSurfaceNormals); 
  // derivs.enroll(mSmoothnessNormalization);             
};

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
// this is from our old implementation where we would just reduce the AV
// based on the velocity and interface normal directions.
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

      const auto ssij = this->pairwiseSurfaceSmoothness(nodeListi,i,nodeListj,j);

      const auto vijhat = (vi-vj).unitVector();
      const auto fi = abs(normi.dot(vijhat));
      const auto fj = abs(normj.dot(vijhat));
      slideCorr = (1.0-ssij) + (ssij)*fi*fj; 
    
    }

    return slideCorr;      

}

//------------------------------------------------------------------------------
// return 1 if pairwise interaction is fully slide and ramp down to zero
// based on the max and min smoothness values for the interaction. These numbers
// were pulled out of a hat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar 
SlideSurface<Dimension>::
pairwiseSurfaceSmoothness(const int nodeListi,
                          const int i, 
                          const int nodeListj,
                          const int j) const {

    typename Dimension::Scalar ssij = 0.0;
    if (this->isSlideSurface(nodeListi,nodeListj)){
      const auto ssMax = max(mSurfaceSmoothness(nodeListi,i),mSurfaceSmoothness(nodeListj,j));
      const auto ssMin = min(mSurfaceSmoothness(nodeListi,i),mSurfaceSmoothness(nodeListj,j));

      //const auto ssijMax = ( 1.0 -  10.0*min(max(0.95-ssMax,0.0),0.10) ); // ramps down 0.98->0.88
      //const auto ssijMin = ( 1.0 -   5.0*min(max(0.85-ssMin,0.0),0.20) ); // ramps down 0.88->0.68
      //ssij = max(ssijMax*ssijMin,0.0);
      //const auto ssijAvg = 0.5*(ssMax+ssMin);
      ssij = ( 1.0 - 10.0*min(max(0.95-(ssMax*ssMin),0.0),0.10) );//ssijMax*ssijMin;//
      //const auto maxFrac = max(mSurfaceFraction(nodeListi,i),mSurfaceFraction(nodeListj,j));
      //ssij = ( 1.0 -  5.0*min(max(maxFrac-0.4,0.0),0.20) );
    }

    return ssij;      

}

//------------------------------------------------------------------------------
// returns pairwise surface normal
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector 
SlideSurface<Dimension>::
pairwiseSurfaceNormal(const int nodeListi,
                      const int i, 
                      const int nodeListj,
                      const int j) const {
    return this->weightedPairwiseSurfaceNormal(nodeListi,i,nodeListj,j,1.0,1.0);      
}


//------------------------------------------------------------------------------
// weighted pairwise smoothness 
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar 
SlideSurface<Dimension>::
weightedPairwiseSurfaceSmoothness(const int nodeListi,
                                  const int i, 
                                  const int nodeListj,
                                  const int j,
                                  const typename Dimension::Scalar weighti,
                                  const typename Dimension::Scalar weightj) const {

    auto ssij = 0.0;

    if (this->isSlideSurface(nodeListi,nodeListj)){

      const auto tiny = std::numeric_limits<double>::epsilon();
      const auto ssi = mSurfaceSmoothness(nodeListi,i);
      const auto ssj = mSurfaceSmoothness(nodeListj,j);

      ssij = (ssi*weighti + ssj * weightj)/max(weighti+weightj,tiny);
      ssij = ( 1.0 - 10.0*min(max(0.95-(ssij),0.0),0.10) );
    }

    return ssij;      

}


//------------------------------------------------------------------------------
// weighted pairwise surface normal 
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector 
SlideSurface<Dimension>::
weightedPairwiseSurfaceNormal(const int nodeListi,
                              const int i, 
                              const int nodeListj,
                              const int j,
                              const typename Dimension::Scalar weighti,
                              const typename Dimension::Scalar weightj) const {

    auto nij = Vector::zero;

    if (this->isSlideSurface(nodeListi,nodeListj)){

      const auto tiny = std::numeric_limits<double>::epsilon();
      const auto smoothnessThreshold = 0.85;

      const auto ni = mSurfaceNormals(nodeListi,i);
      const auto nj = mSurfaceNormals(nodeListj,j);
      const auto ssi = max(mSurfaceSmoothness(nodeListi,i)-smoothnessThreshold,tiny);
      const auto ssj = max(mSurfaceSmoothness(nodeListj,j)-smoothnessThreshold,tiny);

      nij = (ssj*weightj*nj - ssi*weighti*ni);
      nij = nij.unitVector();
    }

    return nij;      

}

//------------------------------------------------------------------------------
// called once on problem start
//------------------------------------------------------------------------------
template<typename Dimension>
void
SlideSurface<Dimension>::
initializeProblemStartup(DataBase<Dimension>& /*dataBase*/){
}

//------------------------------------------------------------------------------
// normals calculated prior to eval derivs step
//------------------------------------------------------------------------------
template<typename Dimension>
void
SlideSurface<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& /*derivs*/,
           typename SlideSurface<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename SlideSurface<Dimension>::ConstBoundaryIterator boundaryEnd,
           const typename Dimension::Scalar /*time*/,
           const typename Dimension::Scalar /*dt*/,
           const TableKernel<Dimension>& W) {

  // if (this->isActive()){

  //   auto  normals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  //   auto  surfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
  //   auto  surfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

  //   surfaceFraction.Zero();
  //   normals.Zero();
  //   surfaceSmoothness.Zero();
    
  //   // computes raw normals (just direction not unit length)
  //   computeSurfaceNormals(W,dataBase,state);

  //   for (ConstBoundaryIterator boundaryItr = boundaryBegin;  
  //        boundaryItr != boundaryEnd; 
  //        ++boundaryItr){ 
  //     (*boundaryItr)->applyFieldListGhostBoundary(normals);
  //     (*boundaryItr)->applyFieldListGhostBoundary(surfaceFraction);
  //     }
  //   for (ConstBoundaryIterator boundaryItr = boundaryBegin; 
  //        boundaryItr != boundaryEnd;
  //        ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  //   //if we are smoothing
  //   if(this->normalsAreSmoothed()){
  //     // apply SPH interpolation to smooth the normal vector field
  //     smoothSurfaceNormals(W,dataBase,state);
  //     for (ConstBoundaryIterator boundaryItr = boundaryBegin;  
  //          boundaryItr != boundaryEnd; 
  //          ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(normals);
  //     for (ConstBoundaryIterator boundaryItr = boundaryBegin; 
  //          boundaryItr != boundaryEnd;
  //          ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  //   }


  //   // get our smoothness metric
  //   computeSurfaceSmoothness(W,dataBase,state);

  //   for (ConstBoundaryIterator boundaryItr = boundaryBegin;
  //        boundaryItr != boundaryEnd;
  //        ++boundaryItr){
  //     (*boundaryItr)->applyFieldListGhostBoundary(surfaceSmoothness);
  //   }
  //   for (ConstBoundaryIterator boundaryItr = boundaryBegin; 
  //        boundaryItr != boundaryEnd;
  //        ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  // } // if active
};  // method



// //------------------------------------------------------------------------------
// // calculate our surface normals
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// SlideSurface<Dimension>::
// computeSurfaceNormals(const TableKernel<Dimension>& W,
//                       const DataBase<Dimension>& dataBase,
//                             State<Dimension>& state) {

//   const auto& connectivityMap = dataBase.connectivityMap();
//   const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
//   const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
//   const auto& massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
//   const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
//         auto  surfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
//         auto  surfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);

//   // how do we define our color field
//   const auto sameMatIsActive = this->surfaceNormalMethod() == SurfaceNormalMethod::SameMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::AllMaterialSurfaceNormals or
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;
//   const auto diffMatIsActive = this->surfaceNormalMethod() == SurfaceNormalMethod::DifferentMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::AllMaterialSurfaceNormals or
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;
//   const auto sameMatCoeff = ( sameMatIsActive ?  1.0 : 0.0);
//   const auto diffMatCoeff = ( diffMatIsActive ? -1.0 : 0.0);
//   const auto isMassWeighted = this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;

//   // Pre-conditions.
//   const auto numNodeLists = massDensity.size();
//   REQUIRE(position.size() == numNodeLists);
//   REQUIRE(mass.size() == numNodeLists);
//   REQUIRE(H.size() == numNodeLists);
//   REQUIRE(surfaceSmoothnessNormalization.size() == numNodeLists);
//   REQUIRE(surfaceNormals.size() == numNodeLists);

//   // The set of interacting node pairs.
//   const auto& pairs = connectivityMap.nodePairList();
//   const auto  npairs = pairs.size();
//   const auto  tiny = 1.0e-25;

//   // do this with a temporary storage variable for now
//   FieldList<Dimension, Tensor> M(FieldStorageType::CopyFields);
//   if(this->gradientsAreCorrected()){
//     for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
//       M.appendNewField("MCorr", massDensity[nodeListi]->nodeList(), Tensor::zero);
//     }
//   }

//   // Now the pair contributions.
// #pragma omp parallel
//   {
//     int i, j, nodeListi, nodeListj;
//     auto surfaceNormals_thread = surfaceNormals.threadCopy();
//     auto surfaceFraction_thread = surfaceFraction.threadCopy();
//     auto M_thread = M.threadCopy();
// #pragma omp for
//     for (auto k = 0u; k < npairs; ++k) {
//       i = pairs[k].i_node;
//       j = pairs[k].j_node;
//       nodeListi = pairs[k].i_list;
//       nodeListj = pairs[k].j_list;

//       const auto sameMatij = (nodeListi == nodeListj);
//       const auto materialCoeff = (sameMatij ? sameMatCoeff : diffMatCoeff);

//       // State for node i
//       const auto& ri = position(nodeListi, i);
//       const auto  mi = mass(nodeListi, i);
//       const auto  rhoi = massDensity(nodeListi, i);
//       const auto& Hi = H(nodeListi, i);
//       const auto  Hdeti = Hi.Determinant();
//       const auto  voli = mi/rhoi;

//       // State for node j
//       const auto& rj = position(nodeListj, j);
//       const auto  mj = mass(nodeListj, j);
//       const auto  rhoj = massDensity(nodeListj, j);
//       const auto& Hj = H(nodeListj, j);
//       const auto  Hdetj = Hj.Determinant();
//       const auto  volj = mj/rhoj;

//       // Kernel weighting and gradient.
//       const auto rij = ri - rj;
//       const auto etai = Hi*rij;
//       const auto etaj = Hj*rij;
//       const auto etaMagi = etai.magnitude();
//       const auto etaMagj = etaj.magnitude();
//       const auto Hetai = Hi*etai.unitVector();
//       const auto Hetaj = Hj*etaj.unitVector();

//       const auto gWi = W.gradValue(etaMagi, Hdeti);
//       const auto gWj = W.gradValue(etaMagj, Hdetj);
//       const auto gradWi = gWi*Hetai;
//       const auto gradWj = gWj*Hetaj;
//       const auto gradWij = (gradWi+gradWj)*0.5;
      
//       const auto Wi = W.kernelValue(etaMagi, Hdeti);
//       const auto Wj = W.kernelValue(etaMagj, Hdetj); 
//       const auto Wij = 0.5*(Wi+Wj);
        
//       if (!sameMatij){
//         surfaceFraction_thread(nodeListi, i) += volj * Wij;
//         surfaceFraction_thread(nodeListj, j) += voli * Wij;
//       }

//       surfaceNormals_thread(nodeListi, i) -=  materialCoeff * (isMassWeighted ? mj : volj) * gradWij;
//       surfaceNormals_thread(nodeListj, j) +=  materialCoeff * (isMassWeighted ? mi : voli) * gradWij;

//       if(this->gradientsAreCorrected()){
//         const auto Mij = rij.dyad(gradWij);
//         M_thread(nodeListi,i) -= abs(materialCoeff) * volj*Mij;
//         M_thread(nodeListj,j) -= abs(materialCoeff) * voli*Mij;
//       }

//     }   // pair loop

// #pragma omp critical
//     {
//       surfaceNormals_thread.threadReduce();
//       surfaceFraction_thread.threadReduce();
//       M_thread.threadReduce();
//     }
//   }   // omp region

//     for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
//      const auto n = surfaceFraction[nodeListi]->numInternalElements();
//  #pragma omp parallel for
//       for (auto i = 0u; i < n; ++i) {

//         //if (surfaceFraction(nodeListi,i)>tiny){
        
//           if(this->gradientsAreCorrected()){
//             const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
//             auto& Mi = M(nodeListi, i);
//             const auto Mdeti = Mi.Determinant();
//             const auto goodM = ( Mdeti > 1.0e-2 and numNeighborsi > Dimension::pownu(2));
//             Mi =  (goodM ? Mi.Inverse() : Tensor::one);
//             surfaceNormals(nodeListi,i) = Mi*surfaceNormals(nodeListi,i);
//           }
//           if(!this->normalsAreSmoothed()){
//             surfaceNormals(nodeListi,i) = (surfaceNormals(nodeListi,i)).unitVector();
//           }
//         //}else{

//         //  surfaceNormals(nodeListi,i) = Vector::zero;

//        // } // if statement
//       }   // node loop
//     }     // nodelist loop
// };        // function



// //------------------------------------------------------------------------------
// // smooth surface normals basis in direction normal to interface
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// SlideSurface<Dimension>::
// smoothSurfaceNormals(const TableKernel<Dimension>& W,
//                      const DataBase<Dimension>& dataBase,
//                            State<Dimension>& state) {

//   const auto& connectivityMap = dataBase.connectivityMap();
//   const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
//   const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
//   const auto& massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
//   const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
//         auto  surfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
//   const auto& surfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);

//   // how do we define our color field
//   const auto sameMatIsActive = this->surfaceNormalMethod() == SurfaceNormalMethod::SameMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::AllMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;
//   const auto diffMatIsActive = this->surfaceNormalMethod() == SurfaceNormalMethod::DifferentMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::AllMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;
//   const auto sameMatCoeff = ( sameMatIsActive ?  1.0 : 0.0);
//   const auto diffMatCoeff = ( diffMatIsActive ? -1.0 : 0.0);
//   const auto isMassWeighted = this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;
  
//   // Pre-conditions.
//   const auto numNodeLists = massDensity.size();
//   REQUIRE(position.size() == numNodeLists);
//   REQUIRE(mass.size() == numNodeLists);
//   REQUIRE(H.size() == numNodeLists);
//   REQUIRE(surfaceSmoothnessNormalization.size() == numNodeLists);
//   REQUIRE(surfaceNormals.size() == numNodeLists);

//   // The set of interacting node pairs.
//   const auto& pairs = connectivityMap.nodePairList();
//   const auto  npairs = pairs.size();
//   const auto  W0 = W.kernelValue(0.0, 1.0);
//   const auto  tiny = 1.0e-25;

//   // do this with a temporary storage variable for now
//   FieldList<Dimension, Vector> n0(FieldStorageType::CopyFields);
//   for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
//     n0.appendNewField("normalsTemp", massDensity[nodeListi]->nodeList(), Vector::zero);
//   }

//   // Now the pair contributions.
// #pragma omp parallel
//   {
//     int i, j, nodeListi, nodeListj;
//     auto surfaceNormals_thread = n0.threadCopy();

// #pragma omp for
//     for (auto k = 0u; k < npairs; ++k) {
//       i = pairs[k].i_node;
//       j = pairs[k].j_node;
//       nodeListi = pairs[k].i_list;
//       nodeListj = pairs[k].j_list;

//       const auto fraci = surfaceFraction(nodeListi,i);
//       const auto fracj = surfaceFraction(nodeListj,j);

//       if(std::min(fraci,fracj)>tiny){
        
//         const auto sameMatij = nodeListi==nodeListj;
//         const auto materialCoeff = (sameMatij ? sameMatCoeff : diffMatCoeff);

//         // State for node i
//         const auto& ri = position(nodeListi, i);
//         const auto  mi = mass(nodeListi, i);
//         const auto  rhoi = massDensity(nodeListi, i);
//         const auto& Hi = H(nodeListi, i);
//         const auto& ni = surfaceNormals(nodeListi,i);

//         // State for node j
//         const auto& rj = position(nodeListj, j);
//         const auto  mj = mass(nodeListj, j);
//         const auto  rhoj = massDensity(nodeListj, j);
//         const auto& Hj = H(nodeListj, j);
//         const auto& nj = surfaceNormals(nodeListj,j);

//         // Kernel weighting and gradient.
//         const auto rij = ri - rj;
//         const auto etai = Hi*rij;
//         const auto etaj = Hj*rij;

//         const auto Wi = W.kernelValue(etai.magnitude(), Hi.Determinant());
//         const auto Wj = W.kernelValue(etaj.magnitude(), Hj.Determinant());

//         surfaceNormals_thread(nodeListi, i) += materialCoeff * (isMassWeighted ? mj : mj/rhoj) * nj * Wi;
//         surfaceNormals_thread(nodeListj, j) += materialCoeff * (isMassWeighted ? mi : mi/rhoi) * ni * Wj;
//       } // if statement
//     }   // pair-loop

// #pragma omp critical
//     {
//       surfaceNormals_thread.threadReduce();
//     }
//   }   // omp reigion

//   // finish with self contribution
//   for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
//      const auto n = surfaceNormals[nodeListi]->numInternalElements();
//  #pragma omp parallel for
//     for (auto i = 0u; i < n; ++i) {
//       if (surfaceFraction(nodeListi,i)>tiny){
//         const auto Hdeti = H(nodeListi,i).Determinant();
//         const auto ni = mass(nodeListi,i)/massDensity(nodeListi,i) * surfaceNormals(nodeListi,i) * Hdeti * W0;
//         surfaceNormals(nodeListi,i) = (ni + n0(nodeListi,i)).unitVector();
//       }else{
//         surfaceNormals(nodeListi,i)=Vector::zero;
//       }
//     }   // node loop
//   }     // nodelist loop
// };      // function


// //------------------------------------------------------------------------------
// // determine the smoothness of our surface
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void
// SlideSurface<Dimension>::
// computeSurfaceSmoothness(const TableKernel<Dimension>& W,
//                          const DataBase<Dimension>& dataBase,
//                                State<Dimension>& state) {

//   const auto& connectivityMap = dataBase.connectivityMap();
//   const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
//   const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
//   const auto& massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
//   const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
//   const auto& surfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
//   const auto& surfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
//         auto  surfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

//   const auto tiny = 1e-25;

//   // how do we define our color field
//   const auto sameMatIsActive = this->surfaceNormalMethod() == SurfaceNormalMethod::SameMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::AllMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;
//   const auto diffMatIsActive = this->surfaceNormalMethod() == SurfaceNormalMethod::DifferentMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::AllMaterialSurfaceNormals or 
//                                this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;
//   const auto sameMatCoeff = ( sameMatIsActive ?  1.0 : 0.0);
//   const auto diffMatCoeff = ( diffMatIsActive ? -1.0 : 0.0);
//   const auto isMassWeighted = this->surfaceNormalMethod() == SurfaceNormalMethod::MassWeightedSurfaceNormals;

//   // Pre-conditions.
//   const auto numNodeLists = massDensity.size();
//   REQUIRE(position.size() == numNodeLists);
//   REQUIRE(mass.size() == numNodeLists);
//   REQUIRE(H.size() == numNodeLists);
//   REQUIRE(surfaceNormals.size() == numNodeLists);
//   REQUIRE(surfaceFraction.size() == numNodeLists);
//   REQUIRE(surfaceSmoothness.size() == numNodeLists);

//   // The set of interacting node pairs.
//   const auto& pairs = connectivityMap.nodePairList();
//   const auto  npairs = pairs.size();

//     // do this with a temporary storage variable for now
//   FieldList<Dimension, Scalar> surfaceSmoothnessNormalization(FieldStorageType::CopyFields);
//   for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
//     surfaceSmoothnessNormalization.appendNewField("smoothnessNormalization", massDensity[nodeListi]->nodeList(), 0.0);
//   }

//   // Now the pair contributions.
// #pragma omp parallel
//   {
//     int i, j, nodeListi, nodeListj;
//     auto surfaceSmoothness_thread = surfaceSmoothness.threadCopy();
//     auto surfaceSmoothnessNormalization_thread = surfaceSmoothnessNormalization.threadCopy();

// #pragma omp for
//     for (auto k = 0u; k < npairs; ++k) {
//       i = pairs[k].i_node;
//       j = pairs[k].j_node;
//       nodeListi = pairs[k].i_list;
//       nodeListj = pairs[k].j_list;

//       //const auto fraci = surfaceFraction(nodeListi,i);
//       //const auto fracj = surfaceFraction(nodeListj,j);

//       const auto sameMatij = (nodeListi == nodeListj);

//       //if(diffMatij){
        
//         const auto materialCoeff = (sameMatij ? sameMatCoeff : diffMatCoeff);

//         // State for node i
//         const auto& ri = position(nodeListi, i);
//         const auto  mi = mass(nodeListi, i);
//         const auto  rhoi = massDensity(nodeListi, i);
//         const auto& Hi = H(nodeListi, i);
//         const auto& ni = surfaceNormals(nodeListi,i);
//         const auto  voli = mi/rhoi;

//         // State for node j
//         const auto& rj = position(nodeListj, j);
//         const auto  mj = mass(nodeListj, j);
//         const auto  rhoj = massDensity(nodeListj, j);
//         const auto& Hj = H(nodeListj, j);
//         const auto& nj = surfaceNormals(nodeListj,j);
//         const auto  volj = mj/rhoj;

//         // Kernel weighting and gradient.
//         const auto rij = ri - rj;
//         const auto etai = Hi*rij;
//         const auto etaj = Hj*rij;

//         const auto Wij = 0.5*(W.kernelValue(etai.magnitude(), Hi.Determinant())+
//                               W.kernelValue(etaj.magnitude(), Hj.Determinant()));
//         const auto Wi =  (isMassWeighted ? mj : volj) * Wij;
//         const auto Wj =  (isMassWeighted ? mi : voli) * Wij;

//         surfaceSmoothnessNormalization_thread(nodeListi, i) += Wi;
//         surfaceSmoothnessNormalization_thread(nodeListj, j) += Wj;
       
//           const auto nirij =  ni.dot(rij);
//           const auto njrij = -nj.dot(rij);

//           const auto isSubmerged = (sameMatij ? false : (nirij > 0.0 or njrij > 0.0));
//           if (!isSubmerged){
//             const auto ninja = max(ni.dot(materialCoeff*nj),0.0);
//             surfaceSmoothness_thread(nodeListi, i) += ninja * Wi;
//             surfaceSmoothness_thread(nodeListj, j) += ninja * Wj;
//           }
//       //} // if different mat
//     }   // pair loop

// #pragma omp critical
//     {
//       surfaceSmoothnessNormalization_thread.threadReduce();
//       surfaceSmoothness_thread.threadReduce();
//     }
//   }

//   for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
//      const auto n = surfaceSmoothnessNormalization[nodeListi]->numInternalElements();
//  #pragma omp parallel for
//      for (auto i = 0u; i < n; ++i) {
//        surfaceSmoothness(nodeListi,i) = min(1.0, 
//                                         max(0.0, 
//                                                 surfaceSmoothness(nodeListi,i) / 
//                                                 max(surfaceSmoothnessNormalization(nodeListi,i),tiny)
//                                             ));
//      }
    
//    }

// }; // function

} // spheral namespace