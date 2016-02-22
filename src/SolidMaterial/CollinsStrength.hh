//---------------------------------Spheral++----------------------------------//
// CollinsYieldStrength -- Implements a pressure dependent yield strength
// model appropriate for geological materials.
//
// Since this is solely a yield strength model it takes another StrengthModel
// as an argument to compute the shear modulus.  Perhaps at some point we should
// just split up the ideas of what provides shear modulus and yield strength?
// 
//    See Collins, Melosh, Ivanov, 2004 Appendix, MAPS
//
// Created by JMO, Thu Jan 14 16:40:21 PST 2016
//     Based on python implementation by Megan Syal and Cody Raskin
//----------------------------------------------------------------------------//
#ifndef __Spheral_CollinsStrength_hh__
#define __Spheral_CollinsStrength_hh__

#include "StrengthModel.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class CollinsStrength: public StrengthModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;

  // Constructors, destructor.
  CollinsStrength(const StrengthModel<Dimension>& shearModulusModel,  // Used to compute the shear modulus
                  const double mui,                                   // Coefficient of internal friction
                  const double Y0,                                    // Shear strength at zero pressure
                  const double Ym);                                   // von Mises plastic limit
  virtual ~CollinsStrength();

  // Override the required generic interface.
  virtual void shearModulus(FieldSpace::Field<Dimension, Scalar>& shearModulus,
                            const FieldSpace::Field<Dimension, Scalar>& density,
                            const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                            const FieldSpace::Field<Dimension, Scalar>& pressure) const;

  virtual void yieldStrength(FieldSpace::Field<Dimension, Scalar>& yieldStrength,
                             const FieldSpace::Field<Dimension, Scalar>& density,
                             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                             const FieldSpace::Field<Dimension, Scalar>& pressure,
                             const FieldSpace::Field<Dimension, Scalar>& plasticStrain,
                             const FieldSpace::Field<Dimension, Scalar>& plasticStrainRate) const;

  virtual void soundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
                          const FieldSpace::Field<Dimension, Scalar>& density,
                          const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                          const FieldSpace::Field<Dimension, Scalar>& pressure,
                          const FieldSpace::Field<Dimension, Scalar>& fluidSoundSpeed) const;

  // Access the strength parameters.
  const StrengthModel<Dimension>& shearModulusModel() const;
  double mui() const;
  double Y0() const;
  double Ym() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const StrengthModel<Dimension>& mShearModulusModel;
  double mmui, mY0, mYm;

  // No copying or assignment.
  CollinsStrength(const CollinsStrength&);
  CollinsStrength& operator=(const CollinsStrength&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class CollinsStrength;
  }
}

#endif

