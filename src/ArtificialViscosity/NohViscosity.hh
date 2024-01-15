//---------------------------------Spheral++----------------------------------//
// A version of the Monaghan Gingold pairwise viscosity with a hard wired
// suppression of the viscosity outside the shock.
//
// Created by JMO, Tue Sep  4 21:05:43 PDT 2001
//----------------------------------------------------------------------------//
#ifndef NohViscosity_HH
#define NohViscosity_HH

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {

template<typename Dimension>
class NohViscosity: public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  NohViscosity();
  NohViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~NohViscosity();

  // Initialization method for the beginning of a cycle.
  virtual 
  void initialize(const DataBase<Dimension>& dataBase,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
                  const Scalar time,
                  const Scalar dt,
                  const TableKernel<Dimension>& W);

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
