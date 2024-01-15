//---------------------------------Spheral++----------------------------------//
// GradPressureViscosity.  A viscosity based on pairwise estimates of the
// pressure gradient.
//
// Created by JMO, Fri Dec 13 15:37:42 PST 2002
//----------------------------------------------------------------------------//
#ifndef GradPressureViscosity_HH
#define GradPressureViscosity_HH

#include "ArtificialViscosity.hh"

namespace Spheral {


template<typename Dimension>
class GradPressureViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename vector<ArtificialViscosity<Dimension>*>::iterator iterator;
  typedef typename vector<ArtificialViscosity<Dimension>*>::const_iterator const_iterator;

  // Constructors.
  GradPressureViscosity();
  GradPressureViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~GradPressureViscosity();

  // Initialization method for the beginning of a cycle.
  virtual 
  void initialize(const DataBase<Dimension>& dataBase,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
                  const Scalar time,
                  const Scalar dt,
                  const TableKernel<Dimension>& W);

  // Method to calculate and return the viscous acceleration, work, and pressure,
  // all in one step (efficiency and all).
  virtual void viscousEffects(Vector& acceleration,
                              Scalar& work,
                              Scalar& pressure,
                              const NodeIteratorBase<Dimension>& nodeI,
                              const NodeIteratorBase<Dimension>& nodeJ,
                              const Vector& rij, 
                              const Vector& rijUnit,
                              const Vector& vi, const Vector& vj,
                              const Vector& etai, const Vector& etaj,
                              const Scalar ci, const Scalar cj,
                              const Scalar Pi, const Scalar Pj,
                              const Scalar rhoi, const Scalar rhoj,
                              const Scalar hi, const Scalar hj,
                              const Vector& gradW) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // The interpolation kernel.
  const TableKernel<Dimension>* mKernelPtr;
  FieldList<Dimension, Scalar> mWeight;
};

}

#endif
