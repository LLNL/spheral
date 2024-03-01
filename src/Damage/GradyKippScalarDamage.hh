//---------------------------------Spheral++----------------------------------//
// GradyKippScalarDamage -- Implements a scalar damage model on 
// SolidNodeLists based on the Grady & Kipp crack model.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//
// Created by JMO, Mon Sep 20 21:48:53 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_GradyKippScalarDamage_hh__
#define __Spheral_GradyKippScalarDamage_hh__

#include "ScalarDamageModel.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class SolidNodeList;
  template<typename Dimension> class DataBase;
  template<typename Dimension> class TableKernel;
}

namespace Spheral {

template<typename Dimension>
class GradyKippScalarDamage: 
    public ScalarDamageModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors, destructor.
  GradyKippScalarDamage(SolidNodeList<Dimension>& nodeList,
                        FluidNodeList<Dimension>& damagedNodeList,
                        const double kWeibull,
                        const double mWeibull,
                        const double volume,
                        const TableKernel<Dimension>& kernel,
                        const unsigned seed,
                        const double crackGrowthMultiplier = 0.4,
                        const int minFlawsPerNode = 1,
                        const int minTotalFlaws = 1,
                        const double averageFailure = -1.0);
  virtual ~GradyKippScalarDamage();

  // Important local parameters.
  double kWeibull() const;
  double mWeibull() const;
  unsigned seed() const;

  //**************************************************************************
  // Restart methods.
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //**************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  double mkWeibull;
  double mmWeibull;
  unsigned mSeed;

  // No default constructor, copying or assignment.
  GradyKippScalarDamage();
  GradyKippScalarDamage(const GradyKippScalarDamage&);
  GradyKippScalarDamage& operator=(const GradyKippScalarDamage&);

  using ScalarDamageModel<Dimension>::mFlaws;
};

}

#endif

