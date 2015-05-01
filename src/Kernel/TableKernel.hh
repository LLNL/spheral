//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_TableKernel_hh__
#define __Spheral_TableKernel_hh__

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "Kernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace KernelSpace {

// Forward declarations.
template<typename Dimension> class BSplineKernel;
template<typename Dimension> class W4SplineKernel;
template<typename Dimension> class GaussianKernel;
template<typename Dimension> class SuperGaussianKernel;
template<typename Dimension> class PiGaussianKernel;
template<typename Dimension> class SincKernel;
template<typename Dimension> class NSincPolynomialKernel;
template<typename Dimension> class NBSplineKernel;
template<typename Dimension> class HatKernel;
template<typename Dimension> class QuarticSplineKernel;
template<typename Dimension> class QuinticSplineKernel;

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
              const int numPoints = 1000,
              const double hmult = 1.0);

  // Destructor.
  ~TableKernel();

  // Assignment.
  TableKernel& operator=(const TableKernel& rhs);

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(double etaMagnitude, double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(double etaMagnitude, double Hdet) const;

  // Return the second derivative value for a given normalized distance or position.
  double grad2Value(double etaMagnitude, double Hdet) const;

  // Simultaneously return the kernel value and first derivative.
  std::pair<double, double> kernelAndGradValue(double etaMagnitude, double Hdet) const;

  // Look up the kernel and first derivative for a set.
  void kernelAndGradValues(const std::vector<double>& etaMagnitudes,
                           const std::vector<double>& Hdets,
                           std::vector<double>& kernelValues,
                           std::vector<double>& gradValues) const;

  // Return the equivalent number of nodes per smoothing scale implied by the given
  // sum of kernel values.
  double equivalentNodesPerSmoothingScale(const double Wsum) const;

  // Return the equivalent W sum implied by the given number of nodes per smoothing scale.
  double equivalentWsum(const double nPerh) const;

  // Allow read only access to the node per smoothing scale lookup table.
  const std::vector<double>& nperhValues() const;
  const std::vector<double>& WsumValues() const;

  // Return the number of points being used in the table.
  int numPoints() const;

  // Return the table step size in eta.
  double stepSize() const;
  double stepSizeInv() const;

  // Return the lower bound entry in the table for the given normalized radius.
  int lowerBound(double etaMagnitude) const;

  // Test if the kernel is currently valid.
  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  std::vector<double> mKernelValues;
  std::vector<double> mGradValues;
  std::vector<double> mGrad2Values;
  std::vector<double> mAkernel, mBkernel;
  std::vector<double> mAgrad, mBgrad;
  std::vector<double> mAgrad2, mBgrad2;
  int mNumPoints;
  double mStepSize;

  std::vector<double> mNperhValues;
  std::vector<double> mWsumValues;
  double mMinNperh;
  double mMaxNperh;

  // Initialize the tabular kernel with the given kernels data.
  template<typename KernelType>
  void setTableDataForKernel(const KernelType& kernel, 
                             const int numPoints,
                             const bool gradientAsKernel);

  // Method to initialize the delta kernel values.
  void setParabolicCoeffs();

  // Generic parabolic interpolation.
  double parabolicInterp(const double etaMagnitude,
                         const std::vector<double>& table,
                         const std::vector<double>& a,
                         const std::vector<double>& b) const;

  // Initialize the table relating Wsum to nodes per smoothing scale.
  void setNperhValues(const bool scaleTo1D = false);
};

}
}

#ifndef __GCCXML__
#include "TableKernelInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
}

#endif

