//---------------------------------Spheral++----------------------------------//
// iSALEROCKYieldStrength -- Implements a pressure dependent yield strength
// model appropriate for geological materials.
//
// This is based on the ROCK model in iSALE.  Described in a few places, though
// it seems to keep changing slightly.
// 
//    See Collins, Melosh, Ivanov, 2004 Appendix, MAPS
//    Raducan et al. 2020 (projectile shape effects paper)
//
// Created by JMO, Mon Jul 12 15:05:00 PDT 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_iSALEROCKStrength_hh__
#define __Spheral_iSALEROCKStrength_hh__

#include "StrengthModel.hh"

namespace Spheral {

template<typename Dimension>
class iSALEROCKStrength: public StrengthModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor.
  iSALEROCKStrength(const StrengthModel<Dimension>& shearModulusModel,  // Used to compute the shear modulus
                    const double Yi0,                                   // Intact strength at zero pressure     
                    const double Yiinf,                                 // Intact strength at infinite pressure 
                    const double fi,                                    // Intact internal friction coefficient 
                    const double Yd0,                                   // Damaged strength at zero pressure    
                    const double Ydinf,                                 // Damaged strength at infinite pressure
                    const double fd);                                   // Damaged internal friction coefficient
  virtual ~iSALEROCKStrength();

  // Override the required generic interface.
  virtual bool providesSoundSpeed() const override { return mShearModulusModel.providesSoundSpeed(); }

  virtual void shearModulus(Field<Dimension, Scalar>& shearModulus,
                            const Field<Dimension, Scalar>& density,
                            const Field<Dimension, Scalar>& specificThermalEnergy,
                            const Field<Dimension, Scalar>& pressure,
                            const Field<Dimension, SymTensor>& damage) const override;

  virtual void yieldStrength(Field<Dimension, Scalar>& yieldStrength,
                             const Field<Dimension, Scalar>& density,
                             const Field<Dimension, Scalar>& specificThermalEnergy,
                             const Field<Dimension, Scalar>& pressure,
                             const Field<Dimension, Scalar>& plasticStrain,
                             const Field<Dimension, Scalar>& plasticStrainRate,
                             const Field<Dimension, SymTensor>& damage) const override;

  virtual void soundSpeed(Field<Dimension, Scalar>& soundSpeed,
                          const Field<Dimension, Scalar>& density,
                          const Field<Dimension, Scalar>& specificThermalEnergy,
                          const Field<Dimension, Scalar>& pressure,
                          const Field<Dimension, Scalar>& fluidSoundSpeed,
                          const Field<Dimension, SymTensor>& damage) const override;

  // Access the strength parameters.
  const StrengthModel<Dimension>& shearModulusModel() const;
  double Yi0() const;                                         // Intact strength at zero pressure     
  double Yiinf() const;                                       // Intact strength at infinite pressure 
  double fi() const;                                          // Intact internal friction coefficient 
  double Yd0() const;                                         // Damaged strength at zero pressure    
  double Ydinf() const;                                       // Damaged strength at infinite pressure
  double fd() const;                                          // Damaged internal friction coefficient

private:
  //--------------------------- Private Interface ---------------------------//
  const StrengthModel<Dimension>& mShearModulusModel;
  double mYi0, mYiinf, mfi, mYd0, mYdinf, mfd;

  // No copying or assignment.
  iSALEROCKStrength(const iSALEROCKStrength&);
  iSALEROCKStrength& operator=(const iSALEROCKStrength&);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class iSALEROCKStrength;
}

#endif

