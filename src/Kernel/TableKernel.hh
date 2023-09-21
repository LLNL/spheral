//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_TableKernel_hh__
#define __Spheral_TableKernel_hh__

#include "Kernel.hh"
#include "Utilities/QuadraticInterpolator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class TableKernel: public Kernel<Dimension, TableKernel<Dimension> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  template<typename KernelType>
  TableKernel(const KernelType& kernel,
              const unsigned numPoints = 100u);
  TableKernel(const TableKernel<Dimension>& rhs);

  // Destructor.
  virtual ~TableKernel();

  // Assignment.
  TableKernel& operator=(const TableKernel& rhs);

  // Equivalence
  bool operator==(const TableKernel& rhs) const;

  // Return the kernel weight for a given normalized distance or position.
  Scalar kernelValue(const Scalar etaij, const Scalar Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  Scalar gradValue(const Scalar etaij, const Scalar Hdet) const;

  // Return the second derivative value for a given normalized distance or position.
  Scalar grad2Value(const Scalar etaij, const Scalar Hdet) const;

  // Simultaneously return the kernel value and first derivative.
  void kernelAndGrad(const Vector& etaj, const Vector& etai, const SymTensor& H,
                     Scalar& W,
                     Vector& gradW,
                     Scalar& deltaWsum) const;
  void kernelAndGradValue(const Scalar etaij, const Scalar Hdet,
                          Scalar& W,
                          Scalar& gW) const;

  // Look up the kernel and first derivative for a set.
  void kernelAndGradValues(const std::vector<Scalar>& etaijs,
                           const std::vector<Scalar>& Hdets,
                           std::vector<Scalar>& kernelValues,
                           std::vector<Scalar>& gradValues) const;

  // Return the equivalent number of nodes per smoothing scale implied by the given
  // sum of kernel values.
  Scalar equivalentNodesPerSmoothingScale(const Scalar Wsum) const;

  // Return the equivalent W sum implied by the given number of nodes per smoothing scale.
  Scalar equivalentWsum(const Scalar nPerh) const;

  // Allow read only access to the tabular data.
  const std::vector<Scalar>& nperhValues() const;
  const std::vector<Scalar>& WsumValues() const;

  // Number of points in our lookup data
  size_t numPoints() const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  typedef QuadraticInterpolator InterpolatorType;
  InterpolatorType mInterp, mGradInterp, mGrad2Interp;
  size_t mNumPoints;

  // Data for the nperh lookup algorithm.
  std::vector<Scalar> mNperhValues, mWsumValues;
  Scalar mMinNperh, mMaxNperh;

  // Initialize the table relating Wsum to nodes per smoothing scale.
  void setNperhValues(const bool scaleTo1D = false);
};

}

#include "TableKernelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class TableKernel;
}

#endif

