//---------------------------------Spheral++----------------------------------//
// GHLLC -- HLLC solver w/ gravitational source term
//----------------------------------------------------------------------------//
#ifndef __Spheral_GHLLC_hh__
#define __Spheral_GHLLC_hh__

#include "HLLC.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class WaveSpeedBase;
template<typename Dimension> class LimiterBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class GHLLC : public HLLC<Dimension> {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

public:

  GHLLC(LimiterBase<Dimension>& slopeLimiter,
       WaveSpeedBase<Dimension>& waveSpeedBase,
       const bool linearReconstruction,
       const Vector gravitationalAcceleration);

  ~GHLLC();

  // virtual
  // void interfaceState(const int i,
  //                     const int j,
  //                     const int nodelisti,
  //                     const int nodelistj,
  //                     const Vector& ri,
  //                     const Vector& rj,
  //                     const Scalar& rhoi,   
  //                     const Scalar& rhoj, 
  //                     const Scalar& ci,   
  //                     const Scalar& cj, 
  //                     const Scalar& sigmai,    
  //                     const Scalar& sigmaj,
  //                     const Vector& vi,    
  //                     const Vector& vj,
  //                           Scalar& Pstar,
  //                           Vector& vstar,
  //                           Scalar& rhostari,
  //                           Scalar& rhostarj) const override;

  virtual
  void interfaceState(const int i,
                      const int j,
                      const int nodelisti,
                      const int nodelistj,
                      const Vector& ri,
                      const Vector& rj,
                      const Scalar& rhoi,   
                      const Scalar& rhoj, 
                      const Scalar& ci,   
                      const Scalar& cj, 
                      const Scalar& sigmai,    
                      const Scalar& sigmaj,
                      const Vector& vi,    
                      const Vector& vj,
                      const Vector& DpDxi,
                      const Vector& DpDxj,
                      const Tensor& DvDxi,
                      const Tensor& DvDxj,
                            Scalar& Pstar,
                            Vector& vstar,
                            Scalar& rhostari,
                            Scalar& rhostarj) const override;

  virtual
  void interfaceState(const int i,
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
                            Vector& vstar) const override;

  
  void gravitationalAcceleration(const Vector g);
  Vector gravitationalAcceleration() const;

private:
  Vector mGravitationalAcceleration;
};

}

#include "GHLLCInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class GHLLC;
}

#endif

