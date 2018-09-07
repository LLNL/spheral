//---------------------------------Spheral++----------------------------------//
// GaussianCircularQuadRule
//
// A Gaussian quadrature rule for integrals on the intersection of two 
// circular/spherical node support domains. 
//
// Created by JNJ, Sun Jul 11 12:50:51 PDT 2010
//----------------------------------------------------------------------------//

#ifndef __Spheral_FVPMSpace_GaussianCircularQuadRule__
#define __Spheral_FVPMSpace_GaussianCircularQuadRule__

#include "CircularQuadRule.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

//! \class GaussianCircularQuadRule
//! This class implements a Guassian quadrature rule for integrating 
//! functions on intersections of circular/spherical intersections in space.
template<typename Dimension>
class GaussianCircularQuadRule: public CircularQuadRule<Dimension> {
  // No general case!
};

// 1D specialization -- easy!
template<>
class GaussianCircularQuadRule<Dim<1> >: public CircularQuadRule<Dim<1> > {
public:
  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::SymTensor SymTensor;

  //! Constructs a 1D Gaussian quadrature rule with the given kernel and 
  //! N integration points.
  GaussianCircularQuadRule(const TableKernel<Dim<1> >& W, size_t N);

  // Copy construct and assignment operator are compiler-defined.

  //! Destructor. 
  ~GaussianCircularQuadRule();

  // Overridden methods.

  void setDomains(const Vector& x1, 
                  double r1,
                  const Vector& x2,
                  double r2);

  void getPointsAndWeights(std::vector<Vector>& points,
                           std::vector<double>& weights) const;

private:

  // Number of points.
  size_t mN; 

  // Lists of points and weights.
  std::vector<double> mX, mW;

  // Interval endpoints.
  double mA, mB;
};

// 2D specialization -- harder.
template<>
class GaussianCircularQuadRule<Dim<2> >: public CircularQuadRule<Dim<2> > {
public:
  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;

  //! Constructs a 2D Gaussian quadrature rule with the given kernel.
  GaussianCircularQuadRule(const TableKernel<Dim<2> >& W);

  //! Destructor. 
  ~GaussianCircularQuadRule();

  // Overridden methods.

  void setDomains(const Vector& x1, 
                  double r1,
                  const Vector& x2,
                  double r2);

  void getPointsAndWeights(std::vector<Vector>& points,
                           std::vector<double>& weights) const;

};

// 3D specialization -- yuck.
template<>
class GaussianCircularQuadRule<Dim<3> >: public CircularQuadRule<Dim<3> > {
public:
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;

  //! Constructs a 3D Gaussian quadrature rule with the given kernel.
  GaussianCircularQuadRule(const TableKernel<Dim<3> >& W);

  //! Destructor. 
  ~GaussianCircularQuadRule();

  // Overridden methods.

  void setDomains(const Vector& x1, 
                  double r1,
                  const Vector& x2,
                  double r2);

  void getPointsAndWeights(std::vector<Vector>& points,
                           std::vector<double>& weights) const;

};

}

#endif
