//---------------------------------Spheral++----------------------------------//
// A modified form of the Monaghan & Gingold viscosity, extended to tensor 
// formalism.
//
// Created by J. Michael Owen, Mon Sep  2 14:45:35 PDT 2002
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorMonaghanGingoldViscosity__
#define __Spheral_TensorMonaghanGingoldViscosity__

#include "ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension>
class TensorMonaghanGingoldViscosity: public ArtificialViscosity<Dimension, typename Dimension::Tensor> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors and destuctor
  TensorMonaghanGingoldViscosity(const Scalar Clinear,
                                 const Scalar Cquadratic,
                                 const TableKernel<Dimension>& kernel);
  virtual ~TensorMonaghanGingoldViscosity() = default;

  // No default construction, copying, or assignment
  TensorMonaghanGingoldViscosity() = delete;
  TensorMonaghanGingoldViscosity(const TensorMonaghanGingoldViscosity&) = delete;
  TensorMonaghanGingoldViscosity& operator=(const TensorMonaghanGingoldViscosity&) = delete;

  // We need the velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtficialViscosity (particularly DvDx).
  virtual void QPiij(Tensor& QPiij, Tensor& QPiji,      // result for QPi (Q/rho^2)
                     Scalar& Qij, Scalar& Qji,          // result for viscous pressure
                     const unsigned nodeListi, const unsigned i, 
                     const unsigned nodeListj, const unsigned j,
                     const Vector& xi,
                     const SymTensor& Hi,
                     const Vector& etai,
                     const Vector& vi,
                     const Scalar rhoi,
                     const Scalar csi,
                     const Vector& xj,
                     const SymTensor& Hj,
                     const Vector& etaj,
                     const Vector& vj,
                     const Scalar rhoj,
                     const Scalar csj,
                     const FieldList<Dimension, Scalar>& fCl,
                     const FieldList<Dimension, Scalar>& fCq,
                     const FieldList<Dimension, Tensor>& DvDx) const override;

  // Restart methods.
  virtual std::string label() const override { return "TensorMonaghanGingoldViscosity"; }

protected:
  //--------------------------- Protected Interface ---------------------------//
  using ArtificialViscosity<Dimension, Tensor>::mClinear;
  using ArtificialViscosity<Dimension, Tensor>::mCquadratic;
  using ArtificialViscosity<Dimension, Tensor>::mEpsilon2;
  using ArtificialViscosity<Dimension, Tensor>::mBalsaraShearCorrection;
};

}

#endif
