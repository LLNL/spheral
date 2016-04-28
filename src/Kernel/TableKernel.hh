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

  // // Linearly combine with another kernel.
  // template<typename KernelType>
  // void augment(const KernelType& kernel);

  // Return the kernel weight for a given normalized distance or position.
  double kernelValue(const double etaMagnitude, double Hdet) const;

  // Return the gradient value for a given normalized distance or position.
  double gradValue(const double etaMagnitude, const double Hdet) const;

  // Return the second derivative value for a given normalized distance or position.
  double grad2Value(const double etaMagnitude, const double Hdet) const;

  // Simultaneously return the kernel value and first derivative.
  std::pair<double, double> kernelAndGradValue(const double etaMagnitude, const double Hdet) const;

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

  // Look up the f1 and f2 RZ corrections.
  // Note these methods are only supported for 2D kernels -- other dimensions throw an error.
  double f1(const double etaMagnitude) const;
  double f2(const double etaMagnitude) const;
  double gradf1(const double etaMagnitude) const;
  double gradf2(const double etaMagnitude) const;
  void f1Andf2(const double etaMagnitude,
               double& f1,
               double& f2,
               double& gradf1,
               double& gradf2) const;

  // Allow read only access to the tabular data.
  const std::vector<double>& kernelValues() const;
  const std::vector<double>& gradValues() const;
  const std::vector<double>& grad2Values() const;
  const std::vector<double>& nperhValues() const;
  const std::vector<double>& WsumValues() const;
  const std::vector<double>& f1Values() const;
  const std::vector<double>& f2Values() const;
  const std::vector<double>& gradf1Values() const;
  const std::vector<double>& gradf2Values() const;

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
  // Data for the kernel tabulation.
  std::vector<double> mKernelValues;
  std::vector<double> mGradValues;
  std::vector<double> mGrad2Values;
  std::vector<double> mAkernel, mBkernel, mCkernel;
  std::vector<double> mAgrad, mBgrad, mCgrad;
  std::vector<double> mAgrad2, mBgrad2, mCgrad2;
  int mNumPoints;
  double mStepSize;

  // Data for the nperh lookup algorithm.
  std::vector<double> mNperhValues;
  std::vector<double> mWsumValues;
  double mMinNperh;
  double mMaxNperh;

  // Data for tabulating the RZ f1 and f2 corrections.
  std::vector<double> mf1Values, mAf1, mBf1, mCf1;
  std::vector<double> mf2Values, mAf2, mBf2, mCf2;
  std::vector<double> mGradf1Values, mAgradf1, mBgradf1, mCgradf1;
  std::vector<double> mGradf2Values, mAgradf2, mBgradf2, mCgradf2;

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
                         const std::vector<double>& b,
                         const std::vector<double>& c) const;

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

