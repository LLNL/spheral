//---------------------------------Spheral++----------------------------------//
// A straight implementation of the VonNeuman scalar Q for SPH.
//
// Created by JMO, Sun Jan 14 22:59:51 PST 2001
//----------------------------------------------------------------------------//
#ifndef VonNeumanViscosity_HH
#define VonNeumanViscosity_HH

#include "ArtificialViscosity.hh"

namespace Spheral {
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace FileIOSpace {
      class FileIO;
   }
}

namespace Spheral {
namespace ArtificialViscositySpace {

using Spheral::FieldSpace::FieldList;
using Spheral::FileIOSpace::FileIO;

template<typename Dimension>
class VonNeumanViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef NodeIDIterator<Dimension> IDIterator;

  // Constructors.
  VonNeumanViscosity();
  VonNeumanViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~VonNeumanViscosity();

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
                              const IDIterator& nodeI,
                              const IDIterator& nodeJ,
                              const Vector& rij, 
                              const Vector& vi, const Vector& vj,
                              const Vector& etai, const Vector& etaj,
                              const Scalar ci, const Scalar cj,
                              const Scalar rhoi, const Scalar rhoj,
                              const Vector& gradW) const;

  // Access the viscous internal energy.
  const FieldList<Dimension, Scalar>& viscousInternalEnergy() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual void dumpState(FileIO& file, const string& pathName) const;
  virtual void restoreState(const FileIO& file, const string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  FieldList<Dimension, Scalar> mViscousEnergy;
};

//------------------------------------------------------------------------------
// Return the viscous energy field list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
VonNeumanViscosity<Dimension>::
viscousInternalEnergy() const {
  return mViscousEnergy;
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
