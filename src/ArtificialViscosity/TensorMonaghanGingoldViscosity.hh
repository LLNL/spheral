//---------------------------------Spheral++----------------------------------//
// A modified form of the Monaghan & Gingold viscosity, extended to tensor 
// formalism.
//
// Created by J. Michael Owen, Mon Sep  2 14:45:35 PDT 2002
//----------------------------------------------------------------------------//
#ifndef TensorMonaghanGingoldViscosity_HH
#define TensorMonaghanGingoldViscosity_HH

#include "ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension>
class TensorMonaghanGingoldViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  TensorMonaghanGingoldViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~TensorMonaghanGingoldViscosity();

  // Required method to compute the tensor viscous P/rho^2.
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

  // Restart methods.
  virtual std::string label() const { return "TensorMonaghanGingoldViscosity"; }

private:
  //--------------------------- Private Interface ---------------------------//
  TensorMonaghanGingoldViscosity();
  TensorMonaghanGingoldViscosity(const TensorMonaghanGingoldViscosity&);
  TensorMonaghanGingoldViscosity& operator=(const TensorMonaghanGingoldViscosity&) const;
};

}

#endif
