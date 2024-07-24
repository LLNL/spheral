//---------------------------------Spheral++----------------------------------//
// A version of our tensor viscosity specialized for the SVPHFacetedHydro
// algorithm.
//
// Created by J. Michael Owen, Sat Aug 31 13:31:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef TensorSVPHViscosity_HH
#define TensorSVPHViscosity_HH

#include "ArtificalViscosity/ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension>
class TensorSVPHViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename std::vector<Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  // Constructors.
  TensorSVPHViscosity(Scalar Clinear, Scalar Cquadratic, Scalar fslice);

  // Destructor.
  ~TensorSVPHViscosity();

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

  // Access our internal state.
  Scalar fslice() const;
  void fslice(Scalar x);

  const std::vector<Tensor>& DvDx() const;
  const std::vector<Scalar>& shearCorrection() const;
  const std::vector<Tensor>& Qface() const;

  // Restart methods.
  virtual std::string label() const { return "TensorSVPHViscosity"; }

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mfslice;
  std::vector<Tensor> mDvDx;
  std::vector<Scalar> mShearCorrection;
  std::vector<Tensor> mQface;

  TensorSVPHViscosity();
  TensorSVPHViscosity(const TensorSVPHViscosity&);
  TensorSVPHViscosity& operator=(const TensorSVPHViscosity&) const;
};

}

#endif
