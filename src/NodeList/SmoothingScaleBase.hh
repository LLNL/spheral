//---------------------------------Spheral++----------------------------------//
// SmoothingScaleBase
//
// An abstract base class defining the interface for updating/defining the 
// smoothing scale associated with a FluidNodeList.
//
// Created by JMO, Wed Sep 14 13:27:39 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_SmooothingScaleBase__
#define __Spheral_NodeSpace_SmooothingScaleBase__

#include "Geometry/Dimension.hh"
#include "Mesh/Mesh.hh"

#include <cmath>

namespace Spheral {

template<typename Dimension> class ConnectivityMap;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class TableKernel;
class FileIO;

template<typename Dimension>
class SmoothingScaleBase {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructor.
  SmoothingScaleBase();
  SmoothingScaleBase(const SmoothingScaleBase& rhs);
  SmoothingScaleBase& operator=(const SmoothingScaleBase& rhs);
  virtual ~SmoothingScaleBase();

  // Compute the time derivative and ideal H simultaneously for a Field of H's.
  void newSmoothingScaleAndDerivative(const Field<Dimension, SymTensor>& H,
                                      const Field<Dimension, Vector>& position,
                                      const Field<Dimension, Tensor>& DvDx,
                                      const Field<Dimension, Scalar>& zerothMoment,
                                      const Field<Dimension, Vector>& firstMoment,
                                      const Field<Dimension, SymTensor>& secondMomentEta,
                                      const Field<Dimension, SymTensor>& secondMomentLab,
                                      const ConnectivityMap<Dimension>& connectivityMap,
                                      const TableKernel<Dimension>& W,
                                      const Scalar hmin,
                                      const Scalar hmax,
                                      const Scalar hminratio,
                                      const Scalar nPerh,
                                      Field<Dimension, SymTensor>& DHDt,
                                      Field<Dimension, SymTensor>& Hideal) const;

  // Given the volume and target nperh, compute an effective target hmax
  Scalar hmax(const Scalar Vi, const Scalar nPerh) const;

  //*****************************************************************************
  // Required methods for descendents.
  // Time derivative of the smoothing scale.
  virtual SymTensor
  smoothingScaleDerivative(const SymTensor& H,
                           const Vector& pos,
                           const Tensor& DvDx,
                           const Scalar hmin,
                           const Scalar hmax,
                           const Scalar hminratio,
                           const Scalar nPerh) const = 0;
  
  // Return a new H, with limiting based on the old value.
  virtual SymTensor
  newSmoothingScale(const SymTensor& H,
                    const Vector& pos,
                    const Scalar zerothMoment,
                    const Vector& firstMoment,
                    const SymTensor& secondMomentEta,
                    const SymTensor& secondMomentLab,
                    const TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh,
                    const ConnectivityMap<Dimension>& connectivityMap,
                    const unsigned nodeListi,
                    const unsigned i) const = 0;

  // Determine an "ideal" H for the given moments.
  virtual SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Vector& pos,
                      const Scalar zerothMoment,
                      const Vector& firstMoment,
                      const SymTensor& secondMomentEta,
                      const SymTensor& secondMomentLab,
                      const TableKernel<Dimension>& W,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh,
                      const ConnectivityMap<Dimension>& connectivityMap,
                      const unsigned nodeListi,
                      const unsigned i) const = 0;

  // Compute the new H tensors for a tessellation.
  virtual SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Mesh<Dimension>& mesh,
                      const typename Mesh<Dimension>::Zone& zone,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh) const = 0;
                      

  //*****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#include "SmoothingScaleBaseInline.hh"

#endif
