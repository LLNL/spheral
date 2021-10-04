//========================================================
// HLLC approximate Riemann solver for three wave solution 
//========================================================
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

//========================================================
// Constructor
//========================================================
template<typename Dimension>
HLLC<Dimension>::
HLLC(LimiterBase<Dimension>& slopeLimiter,
     WaveSpeedBase<Dimension>& waveSpeed,
     const bool linearReconstruction,
     const GradientType gradType):
  RiemannSolverBase<Dimension>(slopeLimiter,
                               waveSpeed,
                               linearReconstruction,
                               gradType){

}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
HLLC<Dimension>::
~HLLC(){}


//========================================================
// Interface State scalar
//========================================================
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
                     typename Dimension::Scalar& Pstar,
                     typename Dimension::Vector& vstar,
                     typename Dimension::Scalar& rhostari,
                     typename Dimension::Scalar& rhostarj) const{

  Scalar Si, Sj;

  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto& limiter = this->limiter();
  const auto& waveSpeedObject = this->waveSpeed();

  const auto& DpDx0   = this->DpDx();
  const auto& DvDx0   = this->DvDx();

  const auto rij = ri - rj;
  const auto rhatij = rij.unitVector();

  vstar = 0.5*(vi+vj);
  Pstar = 0.5*(Pi+Pj);

  if (ci > 0.0 or cj > 0.0){


    // default to nodal values
    auto v1i = vi;
    auto v1j = vj;

    auto p1i = Pi;
    auto p1j = Pj;

    // linear reconstruction
    if(this->linearReconstruction()){

      const auto xij = 0.5*(rij);

      // gradients
      const auto DvDxi = DvDx0(nodelisti,i);
      const auto DvDxj = DvDx0(nodelistj,j);
      const auto DpDxi = DpDx0(nodelisti,i);
      const auto DpDxj = DpDx0(nodelistj,j);

      // gradients along line of action
      const auto Dpi = DpDxi.dot(xij);
      const auto Dpj = DpDxj.dot(xij);
      const auto Dvi = DvDxi.dot(xij);
      const auto Dvj = DvDxj.dot(xij);
      const auto Dui = Dvi.dot(rhatij);
      const auto Duj = Dvj.dot(rhatij);
      const auto Dp0 = 0.5*(Pi-Pj);
      const auto Du0 = 0.5*(vi-vj).dot(rhatij);

      if (true){
        const auto oneOverDu0 = (sgn(Du0) / std::max(tiny, abs(Du0)));
        const auto ru0i = Dui*oneOverDu0;
        const auto ru0j = Duj*oneOverDu0;
        //const auto ru0i2 = (ru0i < 0.0 ? ru0i : std::max(ru0i,1.0) );
        //const auto ru0j2 = (ru0j < 0.0 ? ru0j : std::max(ru0j,1.0) );
        const auto phi0ui = limiter.slopeLimiter(ru0i);
        const auto phi0uj = limiter.slopeLimiter(ru0j);
        const auto phiu = std::min(phi0ui,phi0uj);

        const auto oneOverDp0 = (sgn(Dp0) / std::max(tiny, abs(Dp0)));
        const auto rp0i = Dpi*oneOverDp0;
        const auto rp0j = Dpj*oneOverDp0;
        //const auto rp0i2 = (rp0i < 0.0 ? rp0i : std::max(rp0i,1.0) );
        //const auto rp0j2 = (rp0j < 0.0 ? rp0j : std::max(rp0j,1.0) );
        const auto phi0pi = limiter.slopeLimiter(rp0i);
        const auto phi0pj = limiter.slopeLimiter(rp0j);
        const auto phip = std::min(phi0pi,phi0pj);

        const auto phi = std::min(phiu,phip);

        //const auto ratioTotalu = std::min(1.0, 2.0/abs(std::max(tiny,(phi*ru0i + phi*ru0j))) );
        //const auto phiTotalu = limiter.slopeLimiter(ratioTotalu) * ratioTotalu;

        //const auto ratioTotalp = std::min(1.0, 2.0/abs(std::max(tiny,(phi*rp0i + phi*rp0j))) );
        //const auto phiTotalp = limiter.slopeLimiter(ratioTotalp) * ratioTotalp;

        //const auto phiTotal = std::min(phiTotalp,phiTotalu);
        v1i = vi -  phiu * Dvi;
        v1j = vj +  phiu * Dvj;
        p1i = Pi -  phip * Dpi;
        p1j = Pj +  phip * Dpj;
        //v1i = vi -  phi * Dvi;
        // v1j = vj +  phi * Dvj;
        // p1i = Pi -  phi * Dpi;
        // p1j = Pj +  phi * Dpj;

      }else{
        const auto rui = Dui/(sgn(Duj)*std::max(tiny, abs(Duj)));
        const auto ruj = Duj/(sgn(Dui)*std::max(tiny, abs(Dui)));
        const auto xu = std::min(rui,ruj);
        const auto phiu = limiter.slopeLimiter(xu);

        const auto rpi = Dpi/(sgn(Dpj)*std::max(tiny, abs(Dpj)));
        const auto rpj = Dpj/(sgn(Dpi)*std::max(tiny, abs(Dpi)));
        const auto xp = std::min(rpi,rpj);
        const auto phip = limiter.slopeLimiter(xp);

        const auto phi = std::min(phip,phiu);

        v1i = vi - phiu * Dvi;
        v1j = vj + phiu * Dvj;
        p1i = Pi - phip * Dpi;
        p1j = Pj + phip * Dpj;
      } 
  
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

  } // if statement for sound speed check
}// Scalar interface class

template<typename Dimension>
void
HLLC<Dimension>::
interfaceState(const int i,
               const int j,
               const int nodelisti,
               const int nodelistj,
               const Vector& ri,
               const Vector& rj,
               const Scalar& rhoi,   
               const Scalar& rhoj, 
               const Scalar& ci,   
               const Scalar& cj,
               const Scalar& Pi,    
               const Scalar& Pj,
               const Vector& vi,    
               const Vector& vj,
               const SymTensor& Si,    
               const SymTensor& Sj,
               const Tensor& Di,    
               const Tensor& Dj,
                     Vector& Tstar,
                     Vector& vstar) const{



}

} // spheral namespace