//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_TableKernel_hh__
#define __Spheral_TableKernel_hh__

#include "Kernel.hh"
#include "Geometry/Dimension.hh"

#include <vector>
#include <boost/math/interpolators/cardinal_quintic_b_spline.hpp>

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
              const int numPoints = 1000,
              const double hmult = 1.0);

  // Destructor.
  virtual ~TableKernel();

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
  const std::vector<double>& nperhValues() const;
  const std::vector<double>& WsumValues() const;

  // Number of points in our lookup data
  size_t numPoints() const;

  // Test if the kernel is currently valid.
  virtual bool valid() const override;

private:
  //--------------------------- Private Interface ---------------------------//
  // Data for the kernel tabulation.
  typedef boost::math::interpolators::cardinal_quintic_b_spline<double> InterpolatorType;
  InterpolatorType mInterpolator;
  size_t mNumPoints;

  // Data for the nperh lookup algorithm.
  std::vector<double> mNperhValues, mWsumValues;
  double mMinNperh, mMaxNperh;

  // Data for tabulating the RZ f1 and f2 corrections.
  InterpolatorType mf1Interp, mf2Interp;

  // Initialize the tabular kernel with the given kernels data.
  template<typename KernelType>
  void setTableDataForKernel(const KernelType& kernel, 
                             const int numPoints,
                             const bool gradientAsKernel);

  // Method to initialize the delta kernel values.
  void setParabolicCoeffs(const std::vector<double>& table,
                          std::vector<double>& a,
                          std::vector<double>& b,
                          std::vector<double>& c) const;

  // Generic parabolic interpolation.
  double parabolicInterp(const double etaMagnitude,
                         const std::vector<double>& a,
                         const std::vector<double>& b,
                         const std::vector<double>& c) const;

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

