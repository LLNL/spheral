//---------------------------------Spheral++----------------------------------//
// FixedSmoothingScale
//
// Implements the static fixed smoothing scale option.
//
// Created by JMO, Wed Sep 14 13:50:49 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_FixedSmooothingScale__
#define __Spheral_NodeSpace_FixedSmooothingScale__

#include "SmoothingScaleBase.hh"

namespace Spheral {

template<typename Dimension>
class FixedSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructor.
  FixedSmoothingScale();
  FixedSmoothingScale(const FixedSmoothingScale& rhs);
  FixedSmoothingScale& operator=(const FixedSmoothingScale& rhs);
  virtual ~FixedSmoothingScale();

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
                    const unsigned i) const;

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
                      const unsigned i) const;

  // Compute the new H tensors for a tessellation.
  virtual SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Mesh<Dimension>& mesh,
                      const typename Mesh<Dimension>::Zone& zone,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh) const;
};

}

#endif
