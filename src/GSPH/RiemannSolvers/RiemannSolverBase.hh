//---------------------------------Spheral++----------------------------------//
// RiemannSolverBase
//----------------------------------------------------------------------------//
#ifndef __Spheral_RiemannSolverBase_hh__
#define __Spheral_RiemannSolverBase_hh__

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension> class WaveSpeedBase;
template<typename Dimension> class SlopeLimiterBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class RiemannSolverBase {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;

public:

  RiemannSolverBase(SlopeLimiterBase<Dimension>& slopeLimiter,
                    WaveSpeedBase<Dimension>& waveSpeedBase);

  ~RiemannSolverBase();

  virtual 
  void initialize(const DataBase<Dimension>& dataBase,
                  const State<Dimension>& state,
                  const StateDerivatives<Dimension>& derivs,
                  const Scalar time,
                  const Scalar dt,
                  const TableKernel<Dimension>& W);

  virtual
  void interfaceState(const Scalar  Si,   
                      const Scalar  Sj, 
                      const Scalar  sigmai,    
                      const Scalar  sigmaj,
                      const Scalar  ui,    
                      const Scalar  uj,
                            Scalar& Pstar,
                            Scalar& ustar) const = 0;

  const SlopeLimiterBase<Dimension>& slopeLimiter() const;
  void slopeLimiter(SlopeLimiterBase<Dimension>& slopeLimiter);

  const WaveSpeedBase<Dimension>& waveSpeed() const;
  void waveSpeed(WaveSpeedBase<Dimension>& waveSpeed);

  const FieldList<Dimension,Vector>& DpDx() const;
  const FieldList<Dimension,Tensor>& DvDx() const;
private:
  
  SlopeLimiterBase<Dimension>& mSlopeLimiter;   
  WaveSpeedBase<Dimension>& mWaveSpeed;

  FieldList<Dimension, Vector> mDpDx;
  FieldList<Dimension, Tensor> mDvDx;

};

}

#include "RiemannSolverBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class RiemannSolverBase;
}

#endif

