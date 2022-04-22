//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#ifndef LimitedMonaghanGingoldViscosity_HH
#define LimitedMonaghanGingoldViscosity_HH

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {

template<typename Dimension>
class LimitedMonaghanGingoldViscosity: public MonaghanGingoldViscosity<Dimension> {
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
  LimitedMonaghanGingoldViscosity(const Scalar Clinear,
                                  const Scalar Cquadratic,
                                  const bool linearInExpansion,
                                  const bool quadraticInExpansion,
                                  const Scalar etaCritFrac,
                                  const Scalar etaFoldFrac);

  // Destructor.
  virtual ~LimitedMonaghanGingoldViscosity();

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

  // Access the fractions setting the critical spacing for kicking the
  // viscosity back on full force.
  double etaCritFrac() const;
  void etaCritFrac(double val);

  double etaFoldFrac() const;
  void etaFoldFrac(double val);

  // Restart methods.
  virtual std::string label() const { return "LimitedMonaghanGingoldViscosity"; }

protected:
  //--------------------------- Private Interface ---------------------------//
  double mEtaCritFrac, mEtaFoldFrac, mEtaCrit, mEtaFold;
  FieldList<Dimension, Tensor> mGradVel;

private:
  //--------------------------- Private Interface ---------------------------//
  LimitedMonaghanGingoldViscosity();
  LimitedMonaghanGingoldViscosity(const LimitedMonaghanGingoldViscosity&);
  LimitedMonaghanGingoldViscosity& operator=(const LimitedMonaghanGingoldViscosity&) const;
};

}

#else

namespace Spheral {
  // Forward declaration.
  template<typename Dimension> class LimitedMonaghanGingoldViscosity;
}

#endif
