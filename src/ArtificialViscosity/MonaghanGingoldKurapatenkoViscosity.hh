//---------------------------------Spheral++----------------------------------//
// A version of the Monaghan Gingold pairwise viscosity with a hard wired
// suppression of the viscosity outside the shock.
//
// Created by $USER$, $DATE$
//----------------------------------------------------------------------------//
#ifndef MonaghanGingoldKurapatenkoViscosity_HH
#define MonaghanGingoldKurapatenkoViscosity_HH

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {

template<typename Dimension>
class MonaghanGingoldKurapatenkoViscosity: 
    public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  MonaghanGingoldKurapatenkoViscosity();
  MonaghanGingoldKurapatenkoViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~MonaghanGingoldKurapatenkoViscosity();

  // Calculate the artificial internal energy term.
  virtual Scalar viscousInternalEnergy(const NodeIDIterator<Dimension>& nodeI,
                                       const NodeIDIterator<Dimension>& nodeJ,
                                       const Vector& rij, const Vector& vij,
                                       const Vector& etai, const Vector& etaj,
                                       const Scalar ci, const Scalar cj) const;

private:
  //--------------------------- Private Interface ---------------------------//
  FieldList<Dimension, Vector> mNodePositions;
  FieldList<Dimension, SymTensor> mHfield;
  Scalar mCurrentTime;
};

}

#endif
