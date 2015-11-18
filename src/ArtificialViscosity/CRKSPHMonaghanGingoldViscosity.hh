//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is specialized for use with CRKSPH.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#ifndef CRKSPHMonaghanGingoldViscosity_HH
#define CRKSPHMonaghanGingoldViscosity_HH

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

template<typename Dimension>
class CRKSPHMonaghanGingoldViscosity: public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  CRKSPHMonaghanGingoldViscosity(const Scalar Clinear,
                                 const Scalar Cquadratic,
                                 const bool linearInExpansion,
                                 const bool quadraticInExpansion);

  // Destructor.
  virtual ~CRKSPHMonaghanGingoldViscosity();

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
  virtual std::string label() const { return "CRKSPHMonaghanGingoldViscosity"; }

private:
  //--------------------------- Private Interface ---------------------------//
  FieldSpace::FieldList<Dimension, Tensor> mGradVel;

  CRKSPHMonaghanGingoldViscosity();
  CRKSPHMonaghanGingoldViscosity(const CRKSPHMonaghanGingoldViscosity&);
  CRKSPHMonaghanGingoldViscosity& operator=(const CRKSPHMonaghanGingoldViscosity&) const;
};

}
}

#else

namespace Spheral {
  namespace ArtificialViscositySpace {
    // Forward declaration.
    template<typename Dimension> class CRKSPHMonaghanGingoldViscosity;
  }
}

#endif
