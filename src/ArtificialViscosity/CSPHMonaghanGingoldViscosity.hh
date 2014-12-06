//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is specialized for use with CSPH.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#ifndef CSPHMonaghanGingoldViscosity_HH
#define CSPHMonaghanGingoldViscosity_HH

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

template<typename Dimension>
class CSPHMonaghanGingoldViscosity: public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  CSPHMonaghanGingoldViscosity(const Scalar Clinear,
                               const Scalar Cquadratic,
                               const bool linearInExpansion,
                               const bool quadraticInExpansion);

  // Destructor.
  virtual ~CSPHMonaghanGingoldViscosity();

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual void initialize(const DataBaseSpace::DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          ConstBoundaryIterator boundaryBegin,
                          ConstBoundaryIterator boundaryEnd,
                          const Scalar time,
                          const Scalar dt,
                          const KernelSpace::TableKernel<Dimension>& W);

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
  virtual std::string label() const { return "CSPHMonaghanGingoldViscosity"; }

private:
  //--------------------------- Private Interface ---------------------------//
  FieldSpace::FieldList<Dimension, Tensor> mGradVel;

  CSPHMonaghanGingoldViscosity();
  CSPHMonaghanGingoldViscosity(const CSPHMonaghanGingoldViscosity&);
  CSPHMonaghanGingoldViscosity& operator=(const CSPHMonaghanGingoldViscosity&) const;
};

}
}

#else

namespace Spheral {
  namespace ArtificialViscositySpace {
    // Forward declaration.
    template<typename Dimension> class CSPHMonaghanGingoldViscosity;
  }
}

#endif
