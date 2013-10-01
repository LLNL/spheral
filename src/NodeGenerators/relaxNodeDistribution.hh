//------------------------------------------------------------------------------
// Centroidally relax a node distribution in a boundary.
// Optionally the user can specify a weighting function for the nodes.
//------------------------------------------------------------------------------
#ifndef __Spheral_relaxNodeDistribution__
#define __Spheral_relaxNodeDistribution__

#include <vector>

#include "DataBase/DataBase.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

// Base class for weighting the nodes.
template<typename Dimension>
struct WeightingFunctor {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  WeightingFunctor() {}
  virtual ~WeightingFunctor() {}
  virtual void operator()(const Vector& pos,
                          const FacetedVolume& boundary,
                          Vector& posImage,
                          Scalar& weight,
                          Scalar& weightImage) const { 
    this->__call__(pos, boundary, posImage, weight, weightImage);
  }
  virtual void __call__(const Vector& pos,
                        const FacetedVolume& boundary,
                        Vector& posImage,
                        Scalar& weight,
                        Scalar& weightImage) const { 
    const Vector& posb = boundary.closestPoint(pos);
    posImage = posb + (posb - pos);
    weight = 1.0;
    weightImage = 1.0;
  }
};

template<typename Dimension>
void
relaxNodeDistribution(DataBaseSpace::DataBase<Dimension>& dataBase,
                      const typename Dimension::FacetedVolume& boundary,
                      const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                      const WeightingFunctor<Dimension> weighting,
                      const int maxIterations,
                      const double tolerance);
}

#endif
