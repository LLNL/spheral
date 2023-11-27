//---------------------------------Spheral++----------------------------------//
// PiecewiseLinearPorousStrengthModel
//
// Shear modulus and yield are treated as piecewise linear functions of the
// the distension.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_PiecewiseLinearPorousStrengthModel_hh__
#define __Spheral_PiecewiseLinearPorousStrengthModel_hh__

#include "PorousStrengthModel.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <vector>

namespace Spheral {

template<typename Dimension> class StrengthModel;

template<typename Dimension>
class PiecewiseLinearPorousStrengthModel: public PorousStrengthModel<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor.
  PiecewiseLinearPorousStrengthModel(const StrengthModel<Dimension>& solidStrength,
                                     const std::vector<Scalar> porosityAbscissa,
                                     const std::vector<Scalar> shearModulusRatio,
                                     const std::vector<Scalar> yieldStrengthRatio);

  virtual ~PiecewiseLinearPorousStrengthModel();

  // We require the Field only interface!
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

  Scalar getPorousShearModulusRatio(Scalar alpha) const;
  Scalar getPorousYieldStrengthRatio(Scalar alpha) const;
  
  std::vector<Scalar>& porosityAbscissa();
  std::vector<Scalar>& shearModulusRatios();
  std::vector<Scalar>& yieldStrengthRatios();

private:

  std::vector<Scalar> mPorosityAbscissa;
  std::vector<Scalar> mShearModulusRatios;
  std::vector<Scalar> mYieldStrengthRatios;

  // No copying or assignment.
  PiecewiseLinearPorousStrengthModel(const PiecewiseLinearPorousStrengthModel&);
  PiecewiseLinearPorousStrengthModel& operator=(const PiecewiseLinearPorousStrengthModel&);
};

}

#include "PiecewiseLinearPorousStrengthModelInline.hh"

#endif

