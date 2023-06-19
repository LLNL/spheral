//---------------------------------Spheral++----------------------------------//
// SecondOrderArtificialViscosity -- approximate riemann solver
//   Toro E.F., Spruce M., Speares W., (1994) "Restoration of the Contact Surface in
//   the HLL-Riemann Solver," Shock Waves, 4:25-34
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_SecondOrderArtificialViscosity_hh__
#define __Spheral_SecondOrderArtificialViscosity_hh__

#include "RiemannSolverBase.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class WaveSpeedBase;
template<typename Dimension> class LimiterBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class SecondOrderArtificialViscosity : public RiemannSolverBase<Dimension> {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

public:

  SecondOrderArtificialViscosity(const Scalar Cl,
                                 const Scalar Cq,
                                 LimiterBase<Dimension>& slopeLimiter,
                                 WaveSpeedBase<Dimension>& waveSpeedBase,
                                 const bool linearReconstruction);

  ~SecondOrderArtificialViscosity();

  
  virtual
  void interfaceState(const Vector& ri,
                      const Vector& rj,
                      const SymTensor& Hi,
                      const SymTensor& Hj,
                      const Scalar& rhoi,   
                      const Scalar& rhoj, 
                      const Scalar& ci,   
                      const Scalar& cj, 
                      const Scalar& sigmai,    
                      const Scalar& sigmaj,
                      const Vector& vi,    
                      const Vector& vj,
                      const Vector& DrhoDxi,
                      const Vector& DrhoDxj,
                      const Vector& DpDxi,
                      const Vector& DpDxj,
                      const Tensor& DvDxi,
                      const Tensor& DvDxj,
                            Scalar& Pstar,
                            Vector& vstar,
                            Scalar& rhostari,
                            Scalar& rhostarj) const override;


  virtual
  void interfaceState(const Vector& ri,
                      const Vector& rj,
                      const SymTensor& Hi,
                      const SymTensor& Hj,
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
                            Vector& vstar) const override;

  Scalar Cl() const;
  void Cl(const Scalar x);

  Scalar Cq() const;
  void Cq(const Scalar x);

private:
  Scalar mCl;
  Scalar mCq;
  
};

}

#include "SecondOrderArtificialViscosityInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SecondOrderArtificialViscosity;
}

#endif

