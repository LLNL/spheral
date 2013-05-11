//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScale
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 15:01:13 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_ASPHSmooothingScale__
#define __Spheral_NodeSpace_ASPHSmooothingScale__

#include "SmoothingScaleBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace NodeSpace {

template<typename Dimension>
class ASPHSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructor.
  explicit ASPHSmoothingScale();
  ASPHSmoothingScale(const ASPHSmoothingScale& rhs);
  ASPHSmoothingScale& operator=(const ASPHSmoothingScale& rhs);
  virtual ~ASPHSmoothingScale();

  // Time derivative of the smoothing scale.
  virtual
  SymTensor
  smoothingScaleDerivative(const SymTensor& H,
                           const Tensor& DvDx,
                           const Scalar hmin,
                           const Scalar hmax,
                           const Scalar hminratio,
                           const Scalar nPerh,
                           const int maxNumNeighbors) const;
  
  // Return a new H, with limiting based on the old value.
  virtual
  SymTensor
  newSmoothingScale(const SymTensor& H,
                    const Scalar zerothMoment,
                    const SymTensor& secondMoment,
                    const int numNeighbors,
                    const KernelSpace::TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh,
                    const int maxNumNeighbors) const;

  // Determine an "ideal" H for the given moments.
  virtual
  SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Scalar zerothMoment,
                      const SymTensor& secondMoment,
                      const int numNeighbors,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh,
                      const int maxNumNeighbors) const;
};

// We explicitly specialize the time derivatives.
template<> 
Dim<1>::SymTensor
ASPHSmoothingScale<Dim<1> >::smoothingScaleDerivative(const Dim<1>::SymTensor&, 
                                                      const Dim<1>::Tensor&,
                                                      const Dim<1>::Scalar hmin,
                                                      const Dim<1>::Scalar hmax,
                                                      const Dim<1>::Scalar hminratio,
                                                      const Dim<1>::Scalar nPerh,
                                                      const int maxNumNeighbors) const;
template<>
Dim<2>::SymTensor
ASPHSmoothingScale<Dim<2> >::smoothingScaleDerivative(const Dim<2>::SymTensor&, 
                                                      const Dim<2>::Tensor&,
                                                      const Dim<2>::Scalar hmin,
                                                      const Dim<2>::Scalar hmax,
                                                      const Dim<2>::Scalar hminratio,
                                                      const Dim<2>::Scalar nPerh,
                                                      const int maxNumNeighbors) const;
template<>
Dim<3>::SymTensor
ASPHSmoothingScale<Dim<3> >::smoothingScaleDerivative(const Dim<3>::SymTensor&, 
                                                      const Dim<3>::Tensor&,
                                                      const Dim<3>::Scalar hmin,
                                                      const Dim<3>::Scalar hmax,
                                                      const Dim<3>::Scalar hminratio,
                                                      const Dim<3>::Scalar nPerh,
                                                      const int maxNumNeighbors) const;

}
}

#endif
