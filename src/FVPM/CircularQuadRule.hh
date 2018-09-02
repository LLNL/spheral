//---------------------------------Spheral++----------------------------------//
// CircularQuadRule
//
// A base class for computing integrals on the intersection of two 
// circular/spherical node support domains. 
//
// Created by JNJ, Sun Jul 11 12:50:51 PDT 2010
//----------------------------------------------------------------------------//

#ifndef __Spheral_FVPMSpace_CircularQuadRule__
#define __Spheral_FVPMSpace_CircularQuadRule__

#include "QuadRule.hh"

namespace Spheral {

//! \class CircularQuadRule
//! This abstract base defines an interface for quadrature rules approximating
//! integrals on circular/spherical intersections in space.
template<typename Dimension>
class CircularQuadRule: public QuadRule<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  //! Base constructor for the quadrature rule.
  explicit CircularQuadRule(const TableKernel<Dimension>& W);

  // Copy constructor and assignment operator are compiler-defined.

  //! Destructor. 
  ~CircularQuadRule();

  //! Initialize the quadrature rule by defining the two node 
  //! support domains.
  virtual void setDomain(const Vector& x1, 
                         double r1,
                         const Vector& x2,
                         double r2) = 0;

  // Overridden general setDomain method.
  virtual void setDomain(const Vector& x1, 
                         const SymTensor& H1,
                         const Vector& x2,
                         const SymTensor& H2);


};

}

#endif
