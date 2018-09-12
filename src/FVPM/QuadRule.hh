//---------------------------------Spheral++----------------------------------//
// QuadRule
//
// A base class for computing integrals on the intersection of two 
// node support domains.
//
// Created by JNJ, Sun Jul 11 12:50:51 PDT 2010
//----------------------------------------------------------------------------//

#ifndef __Spheral_FVPMSpace_QuadRule__
#define __Spheral_FVPMSpace_QuadRule__

#include <vector>

namespace Spheral {

template<typename Dimension> class TableKernel;

//! \class EllipticIntersectionQuadRule
//! This abstract base defines an interface for quadrature rules approximating
//! integrals on elliptic intersections in space.
template<typename Dimension>
class QuadRule {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  //! Base constructor for the quadrature rule. Needs a kernel to define the 
  //! integration domain in the most general sense.
  explicit QuadRule(const TableKernel<Dimension>& W);

  // Copy constructor and assignment operator are compiler-defined.

  //! Destructor. 
  virtual ~QuadRule();

  //! Sets the integration domain given two nodes and their smoothing tensors.
  virtual void
  setDomain(const Vector& x1,
            const SymTensor& H1,
            const Vector& x2,
            const SymTensor& H2) = 0;

  //! Retrieves the quadrature points and weights that approximate the 
  //! integral on the domain. The given vectors will be resized if necessary.
  //! \param points A vector that will store the quadrature points.
  //! \param weights A vector that will store the quadrature weights.
  virtual void
  getPointsAndWeights(std::vector<Vector>& points,
                      std::vector<double>& weights) const = 0;

protected:
  // Kernel.
  const TableKernel<Dimension>& mW;

private:
  // No default constructor.
  QuadRule();
};

}

#endif
