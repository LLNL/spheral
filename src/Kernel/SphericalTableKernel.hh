//---------------------------------Spheral++----------------------------------//
// SphericalTableKernel
//
// Take a 3D Kernel and build a specialized 1D tabulated version appropriate
// for use with the spherical SPH algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Created by JMO, Wed Dec  2 16:41:20 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalTableKernel_hh__
#define __Spheral_SphericalTableKernel_hh__

#include "TableKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class SphericalTableKernel: public Kernel<Dim<1>, SphericalTableKernel<Dim<1>>> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::Tensor Tensor;
  typedef Dim<1>::SymTensor SymTensor;

  // Constructors.
  SphericalTableKernel(const TableKernel<Dim<3>>& kernel);

  // Destructor.
  virtual ~SphericalTableKernel();

  // Assignment.
  SphericalTableKernel& operator=(const SphericalTableKernel& rhs);

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
  std::vector<double> mAf1, mBf1, mCf1;
  std::vector<double> mAf2, mBf2, mCf2;
  std::vector<double> mAgradf1, mBgradf1, mCgradf1;
  std::vector<double> mAgradf2, mBgradf2, mCgradf2;

  // Initialize the tabular kernel with the given kernels data.
  template<typename KernelType>
  void setSphericalTableDataForKernel(const KernelType& kernel, 
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

#include "SphericalTableKernelInline.hh"

#else

// Forward declaration.
namespace Spheral {
  class SphericalTableKernel;
}

#endif

