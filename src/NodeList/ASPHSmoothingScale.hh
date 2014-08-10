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
                           const Vector& pos,
                           const Tensor& DvDx,
                           const Scalar hmin,
                           const Scalar hmax,
                           const Scalar hminratio,
                           const Scalar nPerh) const;
  
  // Return a new H, with limiting based on the old value.
  virtual
  SymTensor
  newSmoothingScale(const SymTensor& H,
                    const Vector& pos,
                    const Scalar zerothMoment,
                    const SymTensor& secondMoment,
                    const KernelSpace::TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh,
                    const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                    const unsigned nodeListi,
                    const unsigned i) const;

  // Determine an "ideal" H for the given moments.
  virtual
  SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Vector& pos,
                      const Scalar zerothMoment,
                      const SymTensor& secondMoment,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh,
                      const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                      const unsigned nodeListi,
                      const unsigned i) const;

  // Compute the new H tensors for a tessellation.
  virtual SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const MeshSpace::Mesh<Dimension>& mesh,
                      const typename MeshSpace::Mesh<Dimension>::Zone& zone,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh) const;
};

// We explicitly specialize the time derivatives.
template<> 
Dim<1>::SymTensor
ASPHSmoothingScale<Dim<1> >::smoothingScaleDerivative(const Dim<1>::SymTensor&, 
                                                      const Dim<1>::Vector& pos,
                                                      const Dim<1>::Tensor&,
                                                      const Dim<1>::Scalar hmin,
                                                      const Dim<1>::Scalar hmax,
                                                      const Dim<1>::Scalar hminratio,
                                                      const Dim<1>::Scalar nPerh) const;
template<>
Dim<2>::SymTensor
ASPHSmoothingScale<Dim<2> >::smoothingScaleDerivative(const Dim<2>::SymTensor&, 
                                                      const Dim<2>::Vector& pos,
                                                      const Dim<2>::Tensor&,
                                                      const Dim<2>::Scalar hmin,
                                                      const Dim<2>::Scalar hmax,
                                                      const Dim<2>::Scalar hminratio,
                                                      const Dim<2>::Scalar nPerh) const;
template<>
Dim<3>::SymTensor
ASPHSmoothingScale<Dim<3> >::smoothingScaleDerivative(const Dim<3>::SymTensor&, 
                                                      const Dim<3>::Vector& pos,
                                                      const Dim<3>::Tensor&,
                                                      const Dim<3>::Scalar hmin,
                                                      const Dim<3>::Scalar hmax,
                                                      const Dim<3>::Scalar hminratio,
                                                      const Dim<3>::Scalar nPerh) const;

}
}

#endif
