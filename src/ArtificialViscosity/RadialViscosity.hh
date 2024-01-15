//---------------------------------Spheral++----------------------------------//
// A version of the Monaghan Gingold pairwise viscosity with a hard wired
// radial interaction.
//
// Created by $Author, $Date: 2003-07-31 15:44:32 -0700 (Thu, 31 Jul 2003) $
//----------------------------------------------------------------------------//
#ifndef RadialViscosity_HH
#define RadialViscosity_HH

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {

template<typename Dimension>
class RadialViscosity: public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  RadialViscosity();
  RadialViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~RadialViscosity();

  // Calculate the artificial internal energy term.
  virtual Scalar viscousInternalEnergy(const NodeIDIterator<Dimension>& nodeI,
                                       const NodeIDIterator<Dimension>& nodeJ,
                                       const Vector& rij, const Vector& vij,
                                       const Vector& etai, const Vector& etaj,
                                       const Scalar ci, const Scalar cj) const;

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
