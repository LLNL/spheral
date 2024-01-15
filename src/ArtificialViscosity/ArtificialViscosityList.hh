//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- Make a list of Q's which will sum to act like a 
// single Q.  When one viscosity just won't do!
//
// Created by JMO, Thu Jan 18 17:29:33 PST 2001
//----------------------------------------------------------------------------//
#ifndef ArtificialViscosityList_HH
#define ArtificialViscosityList_HH

#include "ArtificialViscosity.hh"

#include <vector>

namespace Spheral {


template<typename Dimension>
class ArtificialViscosityList: public ArtificialViscosity<Dimension>, 
                               public vector<ArtificialViscosity<Dimension>*> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename vector<ArtificialViscosity<Dimension>*>::iterator iterator;
  typedef typename vector<ArtificialViscosity<Dimension>*>::const_iterator const_iterator;

  // Constructors.
  ArtificialViscosityList();
  ArtificialViscosityList(const vector<ArtificialViscosity<Dimension>*>& QPtrs);

  // Destructor.
  ~ArtificialViscosityList();

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
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

  // Add and delete Q's from the list.
  void appendArtificialViscosity(ArtificialViscosity<Dimension>* QPtr);
  void deleteArtificialViscosity(ArtificialViscosity<Dimension>* QPtr);

  // Test if the given Q is in the list.
  bool haveArtificialViscosity(const ArtificialViscosity<Dimension>* QPtr) const;

  // Test if the ArtificialViscosityList is valid.
  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
};
}

#endif
