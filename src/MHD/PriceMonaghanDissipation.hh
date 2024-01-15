//---------------------------------Spheral++----------------------------------//
// A set of dissipation rules for SPMHD used by Price and Monaghan (2004).
//----------------------------------------------------------------------------//
#ifndef PriceMonaghanDissipation_HH
#define PriceMonaghanDissipation_HH

#include "ArtificialViscosity/ArtificialViscosity.hh"

namespace Spheral {

class PriceMonaghanDissipation: public ArtificialViscosity<Dim<3> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  //! Construct a Price/Monaghan artificial dissipation model.
  //! \param alpha The linear artificial viscosity coefficient.
  //! \param alphaU The linear artificial thermal conductivity coefficient.
  //! \param alphaB The linear artificial resistivity coefficient.
  //! \param beta The quadratic coefficient.
  //! \param mu0 The permeability of free space in vacuum.
  PriceMonaghanDissipation(Scalar alpha, 
                           Scalar alphaU,
                           Scalar alphaB,
                           Scalar beta,
                           Scalar mu0);

  //! Destructor.
  ~PriceMonaghanDissipation();

  // Method to apply the dissipative effects.
  void viscousEffects(const DataBase<Dim<3> >& dataBase,
                      const ConnectivityMap<Dim<3> >& connectivityMap,
                      const State<Dim<3> >& state,
                      StateDerivatives<Dim<3> >& derivatives) const;

  // Timestep vote method.
  TimeStepType dt(const DataBase<Dim<3> >& dataBase, 
                  const State<Dim<3> >& state,
                  const StateDerivatives<Dim<3> >& derivs,
                  const Scalar currentTime) const;

  // Test if the thing is valid.
  virtual bool valid() const;

private:
  // No default constructor.
  PriceMonaghanDissipation();

  // Coefficients.
  Scalar mAlpha, mAlphaU, mAlphaB, mBeta, mMu0;

  // Minimum timestep.
  mutable Scalar mMinDt;
};

}

#endif
