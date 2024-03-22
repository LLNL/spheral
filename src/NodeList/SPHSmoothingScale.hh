//---------------------------------Spheral++----------------------------------//
// SPHSmoothingScale
//
// Implements the standard SPH scalar smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 14:55:16 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_SPHSmooothingScale__
#define __Spheral_NodeSpace_SPHSmooothingScale__

#include "SmoothingScaleBase.hh"

#include <float.h>

namespace Spheral {

template<typename Dimension>
class SPHSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor.
  explicit SPHSmoothingScale();
  SPHSmoothingScale(const SPHSmoothingScale& rhs);
  SPHSmoothingScale& operator=(const SPHSmoothingScale& rhs);
  virtual ~SPHSmoothingScale();

  // Time derivative of the smoothing scale.
  virtual
  SymTensor
  smoothingScaleDerivative(const SymTensor& H,
                           const Vector& pos,
                           const Tensor& DvDx,
                           const Scalar hmin,
                           const Scalar hmax,
                           const Scalar hminratio,
                           const Scalar nPerh) const override;
  
  // Return a new H, with limiting based on the old value.
  virtual
  SymTensor
  newSmoothingScale(const SymTensor& H,
                    const FieldList<Dimension, Vector>& pos,
                    const Scalar zerothMoment,
                    const Vector& firstMoment,
                    const TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh,
                    const ConnectivityMap<Dimension>& connectivityMap,
                    const unsigned nodeListi,
                    const unsigned i) const override;

  // Determine an "ideal" H for the given moments.
  virtual
  SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const FieldList<Dimension, Vector>& pos,
                      const Scalar zerothMoment,
                      const Vector& firstMoment,
                      const TableKernel<Dimension>& W,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh,
                      const ConnectivityMap<Dimension>& connectivityMap,
                      const unsigned nodeListi,
                      const unsigned i) const override;

  // Compute the new H tensors for a tessellation.
  virtual SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Mesh<Dimension>& mesh,
                      const typename Mesh<Dimension>::Zone& zone,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh) const override;
};

}

#endif
