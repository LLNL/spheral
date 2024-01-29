//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
// References: 
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// Created by JMO, Sun May 21 23:46:02 PDT 2000
//----------------------------------------------------------------------------//
#ifndef MonaghanGingoldViscosity_HH
#define MonaghanGingoldViscosity_HH

#include "ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension>
class MonaghanGingoldViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  MonaghanGingoldViscosity(const Scalar Clinear,
                           const Scalar Cquadratic,
                           const bool linearInExpansion,
                           const bool quadraticInExpansion);

  // Destructor.
  virtual ~MonaghanGingoldViscosity();

  // The required method to compute the artificial viscous P/rho^2.
  virtual std::pair<Tensor, Tensor> Piij(const unsigned nodeListi, const unsigned i, 
                                         const unsigned nodeListj, const unsigned j,
                                         const Vector& xi,
                                         const Vector& etai,
                                         const Vector& vi,
                                         const Scalar rhoi,
                                         const Scalar csi,
                                         const SymTensor& Hi,
                                         const Vector& xj,
                                         const Vector& etaj,
                                         const Vector& vj,
                                         const Scalar rhoj,
                                         const Scalar csj,
                                         const SymTensor& Hj) const;

  // Access the switches for acting in expansion.
  bool linearInExpansion() const;
  void linearInExpansion(const bool x);

  bool quadraticInExpansion() const;
  void quadraticInExpansion(const bool x);

  // Restart methods.
  virtual std::string label() const { return "MonaghanGingoldViscosity"; }

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mLinearInExpansion, mQuadraticInExpansion;

private:
  //--------------------------- Private Interface ---------------------------//
  MonaghanGingoldViscosity();
  MonaghanGingoldViscosity(const MonaghanGingoldViscosity&);
  MonaghanGingoldViscosity& operator=(const MonaghanGingoldViscosity&) const;
};

}

#endif
