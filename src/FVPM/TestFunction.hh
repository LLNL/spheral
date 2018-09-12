//---------------------------------Spheral++----------------------------------//
// TestFunction
//
// A test function used for the weak form of the fluid equations.
//
// Created by JNJ, Sat Jul 10 12:08:50 PDT 2010
//----------------------------------------------------------------------------//

#ifndef __Spheral_FVPMSpace_TestFunction__
#define __Spheral_FVPMSpace_TestFunction__

#include "QuadRule.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class TableKernel;

//! \class TestFunction
//! This class represents a test function \f$\psi(\vec{x})\f$ used in the
//! discretization of the Finite Volume Point Method (FVPM).
template<typename Dimension>
class TestFunction {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  //! Constructs a test function.
  //! \param W The kernel used to construct this test function.
  //! \param quadRule The quadrature rule used to compute the interaction vectors.
  TestFunction(const TableKernel<Dimension>& W,
               const QuadRule<Dimension>& quadRule);

  //! Destructor. 
  virtual ~TestFunction();

  //! Evaluates the test function \f$\psi\f$, defined in terms of a node \f$i\f$ 
  //! at \f$\vec{x}_i\f$ (with smoothing tensor \f$\mathbf{H}_i\f$) and a set of neighbors 
  //! \f$\{j\}\f$ at positions \f$\{\vec{x}_j\}\f$ (with smoothing tensors \f$\{\mathbf{H}_j\}\f)$.
  //! This function is evaluated at the position \f$\vec{x}\f$.
  //! \param xi The position \f$\vec{x}_i\f$ at which the node \f$i\f$ is located.
  //! \param Hi The smoothing tensor \f$\mathbf{H}_i\f$ defining the spatial influence of the node \f$i\f$.
  //! \param xjs The positions \f$\{\vec{x}_j\}\f$ at which the nodes \f$\{j\}\f$ are located.
  //! \param Hjs The smoothing tensors \f$\{\mathbf{H}_j\}\f$ defining the spatial influences of the nodes \f$\{j\}\f$.
  //! \param x The position \f$\vec{x}\f$ at which the test function \f$\psi\f$ is evaluated.
  //! \returns The value \f$\psi(\vec{x})\f$.
  virtual double
  operator()(const Vector& xi,
             const SymTensor& Hi,
             const std::vector<Vector>& xjs,
             const std::vector<SymTensor>& Hjs,
             const Vector& x) const;

  //! Evaluate the interaction vector 
  //! \f$\vec{\beta}_{ij} \equiv \int_\Omega \frac{W_i(\vec{x}) \nabla W_j(\vec{x}) - W_j(\vec{x}) \nabla W_i\vec{x})}{\left(\sum_k W_k(\vec{x})\right)^2} \text{d}\Omega\f$
  //! for the given nodes \f$i\f$ and \f$j\f$ on the domain \f$\Omega\f$.a
  //! The interaction vector is the approximation of \f$\nabla\psi\f$ integrated 
  //! over the domain.
  //! \param xi The position \f$\vec{x}_i\f$ at which the node \f$i\f$ is located.
  //! \param Hi The smoothing tensor \f$\mathbf{H}_i\f$ defining the spatial influence of the node \f$i\f$.
  //! \param xjs The positions \f$\{\vec{x}_j\}\f$ at which the nodes \f$\{j\}\f$ are located.
  //! \param Hjs The smoothing tensors \f$\{\mathbf{H}_j\}\f$ defining the spatial influences of the nodes \f$\{j\}\f$.
  //! \returns The interaction vector \vec{\beta}_{ij}.
  virtual Vector 
  interactionVector(const Vector& xi,
                    const SymTensor& Hi,
                    const Vector& xj,
                    const SymTensor& Hj,
                    const std::vector<Vector>& xks,
                    const std::vector<SymTensor>& Hks) const;

private:

  // The kernel.
  TableKernel<Dimension>& mKernel;

  // The quadrature rule.
  mutable QuadRule<Dimension> mQuadRule;

  // No copy constructor or assignment operator is allowed.
  TestFunction(const TestFunction&);
  TestFunction& operator=(const TestFunction&);
};

}

#endif
