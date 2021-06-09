//---------------------------------Spheral++----------------------------------//
// A straight implementation of the VonNeuman scalar Q for SPH.
//
// Created by JMO, Sun Jan 14 22:59:51 PST 2001
//----------------------------------------------------------------------------//
#ifndef VonNeumanViscosity_HH
#define VonNeumanViscosity_HH

#include "ArtificialViscosity.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class VonNeumanViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors.
  VonNeumanViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~VonNeumanViscosity();

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual 
  void initialize(const DataBase<Dimension>& dataBase,
                  const State<Dimension>& state,
                  const StateDerivatives<Dimension>& derivs,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
                  typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
                  const Scalar time,
                  const Scalar dt,
                  const TableKernel<Dimension>& W);

  // Require all descendents to return the artificial viscous Pi = P/rho^2 as a tensor.
  // Scalar viscosities should just return a diagonal tensor with their value along the diagonal.
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

  // Access the viscous internal energy.
  const FieldList<Dimension, Scalar>& viscousEnergy() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//
  FieldList<Dimension, Scalar> mViscousEnergy;
};

}

#else

namespace Spheral {
  // Forward declaration.
  template<typename Dimension> class VonNeumanViscosity;
}

#endif
