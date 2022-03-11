//---------------------------------Spheral++----------------------------------//
// HLLC -- approximate riemann solver
//   Toro E.F., Spruce M., Speares W., (1994) "Restoration of the Contact Surface in
//   the HLL-Riemann Solver," Shock Waves, 4:25-34
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

#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
HLLC<Dimension>::
HLLC(LimiterBase<Dimension>& slopeLimiter,
     WaveSpeedBase<Dimension>& waveSpeed,
     const bool linearReconstruction):
  RiemannSolverBase<Dimension>(slopeLimiter,
                               waveSpeed,
                               linearReconstruction){

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
HLLC<Dimension>::
~HLLC(){}


//------------------------------------------------------------------------------
// Interface State scalar
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// HLLC<Dimension>::
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
//                      typename Dimension::Scalar& /*rhostari*/,
//                      typename Dimension::Scalar& /*rhostarj*/) const{

//   Scalar Si, Sj;

//   const auto tiny = std::numeric_limits<Scalar>::epsilon();

//   const auto& limiter = this->limiter();
//   const auto& waveSpeedObject = this->waveSpeed();

//   const auto& DpDx0   = this->DpDx();
//   const auto& DvDx0   = this->DvDx();

//   const auto rij = ri - rj;
//   const auto rhatij = rij.unitVector();

//   vstar = 0.5*(vi+vj);
//   Pstar = 0.5*(Pi+Pj);

//   if (ci > tiny or cj > tiny){


//     // default to nodal values
//     auto v1i = vi;
//     auto v1j = vj;

//     auto p1i = Pi;
//     auto p1j = Pj;

//     // linear reconstruction
//     if(this->linearReconstruction()){

//       // gradients
//       const auto DvDxi = DvDx0(nodelisti,i);
//       const auto DvDxj = DvDx0(nodelistj,j);
//       const auto DpDxi = DpDx0(nodelisti,i);
//       const auto DpDxj = DpDx0(nodelistj,j);

//       // gradients along line of action
//       if (true){
//         this->linearReconstruction(ri,rj, Pi,Pj, DpDxi,DpDxj,
//                                    p1i,p1j);
//         this->linearReconstruction(ri,rj, vi,vj, DvDxi,DvDxj,
//                                    v1i,v1j);
//       }else{

//         const auto xij = 0.5*(rij);   
//         const auto Dpi = DpDxi.dot(xij);
//         const auto Dpj = DpDxj.dot(xij);
//         const auto Dvi = DvDxi.dot(xij);
//         const auto Dvj = DvDxj.dot(xij);
//         const auto Dui = Dvi.dot(rhatij);
//         const auto Duj = Dvj.dot(rhatij);
//         //const auto Dp0 = 0.5*(Pi-Pj);
//         //const auto Du0 = 0.5*(vi-vj).dot(rhatij);
//         const auto rui = Dui/(sgn(Duj)*std::max(tiny, abs(Duj)));
//         const auto ruj = Duj/(sgn(Dui)*std::max(tiny, abs(Dui)));
//         const auto xu = std::min(rui,ruj);
//         const auto phiu = limiter.slopeLimiter(xu);

//         const auto rpi = Dpi/(sgn(Dpj)*std::max(tiny, abs(Dpj)));
//         const auto rpj = Dpj/(sgn(Dpi)*std::max(tiny, abs(Dpi)));
//         const auto xp = std::min(rpi,rpj);
//         const auto phip = limiter.slopeLimiter(xp);

//         v1i = vi - phiu * Dvi;
//         v1j = vj + phiu * Dvj;
//         p1i = Pi - phip * Dpi;
//         p1j = Pj + phip * Dpj;
//       } 
  
//     }

//     const auto ui = v1i.dot(rhatij);
//     const auto uj = v1j.dot(rhatij);
//     const auto wi = v1i - ui*rhatij;
//     const auto wj = v1j - uj*rhatij;

//     waveSpeedObject.waveSpeed(rhoi,rhoj,ci,cj,ui,uj,Si,Sj);

//     const auto denom = safeInv(Si - Sj);

//     const auto ustar = (Si*ui - Sj*uj - p1i + p1j )*denom;
//     const auto wstar = (Si*wi - Sj*wj)*denom;
//     vstar = ustar*rhatij + wstar;
//     Pstar = Sj * (ustar-uj) + p1j;

//   }else{ // if ci & cj too small punt to normal av
//     const auto uij = std::min((vi-vj).dot(rhatij),0.0);
//     Pstar += 0.25 * (rhoi+rhoj) * (uij*uij);
//   }
// }// Scalar interface class


//------------------------------------------------------------------------------
// Interface State fluid hydro
//------------------------------------------------------------------------------
template<typename Dimension>
void
HLLC<Dimension>::
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
                     typename Dimension::Scalar& /*rhostari*/,
                     typename Dimension::Scalar& /*rhostarj*/) const{

  Scalar Si, Sj;

  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto& waveSpeedObject = this->waveSpeed();

  const auto rij = ri - rj;
  const auto rhatij = rij.unitVector();

  vstar = 0.5*(vi+vj);
  Pstar = 0.5*(Pi+Pj);

  if (ci > tiny or cj > tiny){


    // default to nodal values
    auto v1i = vi;
    auto v1j = vj;

    auto p1i = Pi;
    auto p1j = Pj;

    // linear reconstruction
    if(this->linearReconstruction()){

      // gradients along line of action
      this->linearReconstruction(ri,rj, Pi,Pj, DpDxi,DpDxj,
                                 p1i,p1j);
      this->linearReconstruction(ri,rj, vi,vj, DvDxi,DvDxj,
                                 v1i,v1j);
  
    }

    const auto ui = v1i.dot(rhatij);
    const auto uj = v1j.dot(rhatij);
    const auto wi = v1i - ui*rhatij;
    const auto wj = v1j - uj*rhatij;

    waveSpeedObject.waveSpeed(rhoi,rhoj,ci,cj,ui,uj,Si,Sj);

    const auto denom = safeInv(Si - Sj);

    const auto ustar = (Si*ui - Sj*uj - p1i + p1j )*denom;
    const auto wstar = (Si*wi - Sj*wj)*denom;
    vstar = ustar*rhatij + wstar;
    Pstar = Sj * (ustar-uj) + p1j;

  }else{ // if ci & cj too small punt to normal av
    const auto uij = std::min((vi-vj).dot(rhatij),0.0);
    Pstar += 0.25 * (rhoi+rhoj) * (uij*uij);
  }
}// Scalar interface class


template<typename Dimension>
void
HLLC<Dimension>::
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