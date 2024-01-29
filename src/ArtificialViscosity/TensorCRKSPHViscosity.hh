//---------------------------------Spheral++----------------------------------//
// A specialized form of the TensorMonaghanGingoldViscosity for use with CRKSPH.
//
// Created by J. Michael Owen, Wed Nov  5 23:51:31 PST 2014
//----------------------------------------------------------------------------//
#ifndef TensorCRKSPHViscosity_HH
#define TensorCRKSPHViscosity_HH

#include "TensorMonaghanGingoldViscosity.hh"

namespace Spheral {

template<typename Dimension>
class TensorCRKSPHViscosity: public TensorMonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  TensorCRKSPHViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~TensorCRKSPHViscosity();

  // Required method to compute the tensor viscous P/rho^2.
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
  virtual std::string label() const { return "TensorCRKSPHViscosity"; }

protected:
  //--------------------------- Protected Interface ---------------------------//
  virtual void calculateSigmaAndGradDivV(const DataBase<Dimension>& dataBase,
                                         const State<Dimension>& state,
                                         const StateDerivatives<Dimension>& derivs,
                                         const TableKernel<Dimension>& W,
                                         ConstBoundaryIterator boundaryBegin,
                                         ConstBoundaryIterator boundaryEnd);

private:
  //--------------------------- Private Interface ---------------------------//
  TensorCRKSPHViscosity();
  TensorCRKSPHViscosity(const TensorCRKSPHViscosity&);
  TensorCRKSPHViscosity& operator=(const TensorCRKSPHViscosity&) const;

  FieldList<Dimension, Tensor> mGradVel;
};

}

#endif
