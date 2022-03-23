//---------------------------------Spheral++----------------------------------//
// GHLLC -- HLLC with gravitational source term
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"

#include "Hydro/HydroFieldNames.hh"
#include "GSPH/GSPHFieldNames.hh"

#include "GSPH/WaveSpeeds/WaveSpeedBase.hh"
#include "GSPH/Limiters/LimiterBase.hh"
#include "GSPH/RiemannSolvers/HLLC.hh"
#include "GSPH/RiemannSolvers/GHLLC.hh"

#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
GHLLC<Dimension>::
GHLLC(LimiterBase<Dimension>& slopeLimiter,
     WaveSpeedBase<Dimension>& waveSpeed,
     const bool linearReconstruction,
     const typename Dimension::Vector gravitationalAcceleration):
  HLLC<Dimension>(slopeLimiter,
                  waveSpeed,
                  linearReconstruction),
  mGravitationalAcceleration(gravitationalAcceleration){

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GHLLC<Dimension>::
~GHLLC(){}


//------------------------------------------------------------------------------
// Interface State fluid hydro
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// GHLLC<Dimension>::
// interfaceState(const int i,
//                const int j,
//                const int nodelisti,
//                const int nodelistj,
//                const typename Dimension::Vector& ri,
//                const typename Dimension::Vector& rj,
//                const typename Dimension::Scalar& rhoi,   
//                const typename Dimension::Scalar& rhoj, 
//                const typename Dimension::Scalar& ci,   
//                const typename Dimension::Scalar& cj,
//                const typename Dimension::Scalar& Pi,    
//                const typename Dimension::Scalar& Pj,
//                const typename Dimension::Vector& vi,    
//                const typename Dimension::Vector& vj,
//                      typename Dimension::Scalar& Pstar,
//                      typename Dimension::Vector& vstar,
//                      typename Dimension::Scalar& rhostari,
//                      typename Dimension::Scalar& rhostarj) const{

//     // pressure + linear grav contribution
//     const auto rhogh = 0.5*(rhoi+rhoj)*mGravitationalAcceleration.dot(ri-rj);
//     const auto p1i = Pi - rhogh;
//     const auto p1j = Pj + rhogh;
//     HLLC<Dimension>::interfaceState(i,
//                                     j,
//                                     nodelisti, 
//                                     nodelistj,
//                                     ri,
//                                     rj,
//                                     rhoi, 
//                                     rhoj,
//                                     ci,   
//                                     cj,
//                                     p1i,
//                                     p1j,
//                                     vi,
//                                     vj,
//                                     Pstar,
//                                     vstar,
//                                     rhostari,
//                                     rhostarj);
// }// Scalar interface class


template<typename Dimension>
void
GHLLC<Dimension>::
interfaceState(const int i,
               const int j,
               const int nodelisti,
               const int nodelistj,
               const typename Dimension::Vector& ri,
               const typename Dimension::Vector& rj,
               const typename Dimension::Scalar& rhoi,   
               const typename Dimension::Scalar& rhoj, 
               const typename Dimension::Scalar& ci,   
               const typename Dimension::Scalar& cj,
               const typename Dimension::Scalar& Pi,    
               const typename Dimension::Scalar& Pj,
               const typename Dimension::Vector& vi,    
               const typename Dimension::Vector& vj,
               const typename Dimension::Vector& DpDxi,    
               const typename Dimension::Vector& DpDxj,
               const typename Dimension::Tensor& DvDxi,    
               const typename Dimension::Tensor& DvDxj,
                     typename Dimension::Scalar& Pstar,
                     typename Dimension::Vector& vstar,
                     typename Dimension::Scalar& rhostari,
                     typename Dimension::Scalar& rhostarj) const{

    // pressure + linear grav contribution
    const auto rhogh = 0.5*(rhoi+rhoj)*mGravitationalAcceleration.dot(ri-rj);
    const auto p1i = Pi - rhogh;
    const auto p1j = Pj + rhogh;
    HLLC<Dimension>::interfaceState(i,
                                    j,
                                    nodelisti, 
                                    nodelistj,
                                    ri,
                                    rj,
                                    rhoi, 
                                    rhoj,
                                    ci,   
                                    cj,
                                    p1i,
                                    p1j,
                                    vi,
                                    vj,
                                    DpDxi,
                                    DpDxj,
                                    DvDxi,
                                    DvDxj,
                                    Pstar,
                                    vstar,
                                    rhostari,
                                    rhostarj);
}// Scalar interface class

template<typename Dimension>
void
GHLLC<Dimension>::
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