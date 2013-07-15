//---------------------------------Spheral++----------------------------------//
// The plain old von Neuman/Richtmyer viscosity, using a MASH formalism to
// evaluate the divergence of the velocity.
//
// Created by JMO, Wed Sep 21 20:57:00 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_MASHVonNeumanViscosity__
#define __Spheral_MASHVonNeumanViscosity__

#include "ArtificialViscosity.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "FileIO/FileIO.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

template<typename Dimension>
class MASHVonNeumanViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef NodeIteratorBase<Dimension> IDIterator;

  // Constructors.
  MASHVonNeumanViscosity(const Scalar Clinear,
                         const Scalar Cquadratic);

  // Destructor.
  ~MASHVonNeumanViscosity();

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual 
  void initialize(const DataBaseSpace::DataBase<Dimension>& dataBase,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
		  const Scalar time,
		  const Scalar dt,
                  const KernelSpace::TableKernel<Dimension>& W);

  // Method to calculate and return the viscous acceleration, work, and pressure,
  // all in one step (efficiency and all).
  virtual void viscousEffects(Vector& acceleration,
                              Scalar& work,
                              Scalar& pressure,
                              const IDIterator& nodeI,
                              const IDIterator& nodeJ,
                              const Vector& rij, 
                              const Vector& rijUnit,
                              const Vector& vi, const Vector& vj,
                              const Vector& etai, const Vector& etaj,
                              const Scalar ci, const Scalar cj,
                              const Scalar Pi, const Scalar Pj,
                              const Scalar rhoi, const Scalar rhoj,
                              const Scalar hi, const Scalar hj,
                              const Vector& gradW) const;

  // Access the computed value for the divergence of the velocity.
  const FieldSpace::FieldList<Dimension, Scalar>& velocityDivergence() const;

  // Access the MASH gradient correction.
  const FieldSpace::FieldList<Dimension, Tensor>& correction() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  FieldSpace::FieldList<Dimension, Scalar> mVelocityDivergence;
  FieldSpace::FieldList<Dimension, Tensor> mCorrection;

  MASHVonNeumanViscosity();
  MASHVonNeumanViscosity& operator=(const MASHVonNeumanViscosity& rhs);
  
};

//------------------------------------------------------------------------------
// Return the velocity divergence field list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>&
MASHVonNeumanViscosity<Dimension>::
velocityDivergence() const {
  return mVelocityDivergence;
}

//------------------------------------------------------------------------------
// The MASH gradient correction term.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>&
MASHVonNeumanViscosity<Dimension>::
correction() const {
  return mCorrection;
}

}
}

#else

namespace Spheral {
namespace ArtificialViscositySpace {
// Forward declaration.
template<typename Dimension> class VonNeumanViscosity;
}
}

#endif
