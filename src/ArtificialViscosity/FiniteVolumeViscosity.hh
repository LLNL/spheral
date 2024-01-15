//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructred the 
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#ifndef FiniteVolumeViscosity_HH
#define FiniteVolumeViscosity_HH

#include "ArtificialViscosity.hh"
#include "Field/FieldList.hh"

namespace Spheral {

template<typename Dimension>
class FiniteVolumeViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  FiniteVolumeViscosity(const Scalar Clinear,
                        const Scalar Cquadratic,
                        const bool scalar);

  // Destructor.
  virtual ~FiniteVolumeViscosity();

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual void initialize(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          ConstBoundaryIterator boundaryBegin,
                          ConstBoundaryIterator boundaryEnd,
                          const Scalar time, 
                          const Scalar dt,
                          const TableKernel<Dimension>& W);

  // The required method to compute the artificial viscous P/rho^2.
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
  virtual std::string label() const { return "FiniteVolumeViscosity"; }

  // Access the internal state.
  bool scalar() const;
  const FieldList<Dimension, Tensor>& DvDx() const;

private:
  //--------------------------- Private Interface ---------------------------//
  bool mScalar;
  FieldList<Dimension, Tensor> mDvDx;

  FiniteVolumeViscosity();
  FiniteVolumeViscosity(const FiniteVolumeViscosity&);
  FiniteVolumeViscosity& operator=(const FiniteVolumeViscosity&) const;
};

}

#endif
