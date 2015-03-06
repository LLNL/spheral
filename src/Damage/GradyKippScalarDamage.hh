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
  namespace NodeSpace {
    template<typename Dimension> class SolidNodeList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
}

namespace Spheral {
namespace PhysicsSpace {

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
  GradyKippScalarDamage(SolidMaterial::SolidNodeList<Dimension>& nodeList,
                        NodeSpace::FluidNodeList<Dimension>& damagedNodeList,
                        const double kWeibull,
                        const double mWeibull,
                        const double volume,
                        const KernelSpace::TableKernel<Dimension>& kernel,
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
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
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

#ifndef __GCCXML__
  using ScalarDamageModel<Dimension>::mFlaws;
#endif

};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class GradyKippScalarDamage;
  }
}

#endif

