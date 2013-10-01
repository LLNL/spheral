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

namespace Spheral {

// Base class for weighting the nodes.
template<typename Dimension>
struct WeightingFunctor {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  FacetedVolume& boundary;

  WeightingFunctor(const FacetedVolume& b): boundary(b) {}
  virtual ~WeightingFunctor() {}
  virtual double operator()(const Vector& position) { return 1.0; }
};

template<typename Dimension>
void
relaxNodeDistribution(DataBaseSpace::DataBase<Dimension>& dataBase,
                      const typename Dimension::FacetedVolume& boundary,
                      const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                      const WeightingFunctor<Dimension>& weighting,
                      const int maxIterations = 100,
                      const double tolerance = 1.0e-10);
}

#endif
