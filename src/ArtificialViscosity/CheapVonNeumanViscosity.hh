//---------------------------------Spheral++----------------------------------//
// A specialized version of the scalar VonNeuman viscosity, which uses the 
// last stored version of DvDx to compute the velocity divergence.
//
// Created by JMO, Sun Jan 14 22:59:51 PST 2001
//----------------------------------------------------------------------------//
#ifndef CheapVonNeumanViscosity_HH
#define CheapVonNeumanViscosity_HH

#include "ArtificialViscosity.hh"
#include "VonNeumanViscosity.hh"

namespace Spheral {

template<typename Dimension>
class CheapVonNeumanViscosity: public VonNeumanViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  CheapVonNeumanViscosity();
  CheapVonNeumanViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~CheapVonNeumanViscosity();

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual 
  void initialize(const DataBase<Dimension>& dataBase,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
		  const Scalar time,
		  const Scalar dt,
                  const TableKernel<Dimension>& W);

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
