//---------------------------------Spheral++----------------------------------//
// A specialized form of the TensorMonaghanGingoldViscosity for use with CSPH.
//
// Created by J. Michael Owen, Wed Nov  5 23:51:31 PST 2014
//----------------------------------------------------------------------------//
#ifndef TensorCSPHViscosity_HH
#define TensorCSPHViscosity_HH

#include "TensorMonaghanGingoldViscosity.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

template<typename Dimension>
class TensorCSPHViscosity: public TensorMonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  TensorCSPHViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~TensorCSPHViscosity();

  // Restart methods.
  virtual std::string label() const { return "TensorCSPHViscosity"; }

protected:
  //--------------------------- Protected Interface ---------------------------//
  virtual void calculateSigmaAndGradDivV(const DataBaseSpace::DataBase<Dimension>& dataBase,
                                         const State<Dimension>& state,
                                         const StateDerivatives<Dimension>& derivs,
                                         const KernelSpace::TableKernel<Dimension>& W,
                                         ConstBoundaryIterator boundaryBegin,
                                         ConstBoundaryIterator boundaryEnd);

private:
  //--------------------------- Private Interface ---------------------------//
  TensorCSPHViscosity();
  TensorCSPHViscosity(const TensorCSPHViscosity&);
  TensorCSPHViscosity& operator=(const TensorCSPHViscosity&) const;

  FieldSpace::FieldList<Dimension, Tensor> mGradVel;
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace ArtificialViscositySpace {
    template<typename Dimension> class TensorCSPHViscosity;
  }
}

#endif
