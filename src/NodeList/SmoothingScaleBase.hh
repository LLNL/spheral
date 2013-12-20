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

#include "Mesh/Mesh.hh"

namespace Spheral {

namespace NeighborSpace {
  template<typename Dimension> class ConnectivityMap;
}
namespace FieldSpace {
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;
}
namespace KernelSpace {
  template<typename Dimension> class TableKernel;
}
namespace FileIOSpace {
  class FileIO;
}

namespace NodeSpace {

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
  void newSmoothingScaleAndDerivative(const FieldSpace::Field<Dimension, SymTensor>& H,
                                      const FieldSpace::Field<Dimension, Tensor>& DvDx,
                                      const FieldSpace::Field<Dimension, Scalar>& zerothMoment,
                                      const FieldSpace::Field<Dimension, SymTensor>& secondMoment,
                                      const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                                      const KernelSpace::TableKernel<Dimension>& W,
                                      const Scalar hmin,
                                      const Scalar hmax,
                                      const Scalar hminratio,
                                      const Scalar nPerh,
                                      const int maxNumNeighbors,
                                      FieldSpace::Field<Dimension, SymTensor>& DHDt,
                                      FieldSpace::Field<Dimension, SymTensor>& Hideal) const;

  //*****************************************************************************
  // Required methods for descendents.
  // Time derivative of the smoothing scale.
  virtual SymTensor
  smoothingScaleDerivative(const SymTensor& H,
                           const Tensor& DvDx,
                           const Scalar hmin,
                           const Scalar hmax,
                           const Scalar hminratio,
                           const Scalar nPerh,
                           const int maxNumNeighbors) const = 0;
  
  // Return a new H, with limiting based on the old value.
  virtual SymTensor
  newSmoothingScale(const SymTensor& H,
                    const Scalar zerothMoment,
                    const SymTensor& secondMoment,
                    const int numNeighbors,
                    const KernelSpace::TableKernel<Dimension>& W,
                    const Scalar hmin,
                    const Scalar hmax,
                    const Scalar hminratio,
                    const Scalar nPerh,
                    const int maxNumNeighbors) const = 0;

  // Determine an "ideal" H for the given moments.
  virtual SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const Scalar zerothMoment,
                      const SymTensor& secondMoment,
                      const int numNeighbors,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const Scalar hmin,
                      const Scalar hmax,
                      const Scalar hminratio,
                      const Scalar nPerh,
                      const int maxNumNeighbors) const = 0;

  // Compute the new H tensors for a tessellation.
  virtual SymTensor
  idealSmoothingScale(const SymTensor& H,
                      const MeshSpace::Mesh<Dimension>& mesh,
                      const typename MeshSpace::Mesh<Dimension>::Zone& zone,
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
}

#endif
