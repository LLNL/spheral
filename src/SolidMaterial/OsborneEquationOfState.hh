//---------------------------------Spheral++----------------------------------//
// OsborneEquationOfState -- Osborne  equation of state.
//
// Reference: PAGOSA Physics manual, LA-14425-M
//----------------------------------------------------------------------------//
#ifndef OsborneEquationOfState_HH
#define OsborneEquationOfState_HH

#include <limits>
#include "SolidEquationOfState.hh"

// Forward declarations.
namespace Spheral {
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }
}

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class OsborneEquationOfState: public SolidEquationOfState<Dimension> { 

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  OsborneEquationOfState(const double referenceDensity,
                         const double etamin,
                         const double etamax,
                         const double a1,
                         const double a2pos,
                         const double a2neg,
                         const double b0,
                         const double b1,
                         const double b2pos,
                         const double b2neg,
                         const double c0,
                         const double c1,
                         const double c2pos,
                         const double c2neg,
                         const double E0,
                         const double atomicWeight,
                         const Material::PhysicalConstants& constants,
                         const double externalPressure,
                         const double minimumPressure,
                         const double maximumPressure,
                         const Material::MaterialPressureMinType minPressureType);
  virtual ~OsborneEquationOfState();

  // We require any equation of state to define the following methods for Fields.
  virtual void setPressure(FieldSpace::Field<Dimension, Scalar>& pressure,
                           const FieldSpace::Field<Dimension, Scalar>& massDensity,
                           const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setTemperature(FieldSpace::Field<Dimension, Scalar>& temperature,
                              const FieldSpace::Field<Dimension, Scalar>& massDensity,
                              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setSpecificThermalEnergy(FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                                        const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                        const FieldSpace::Field<Dimension, Scalar>& temperature) const;

  virtual void setSpecificHeat(FieldSpace::Field<Dimension, Scalar>& specificHeat,
                               const FieldSpace::Field<Dimension, Scalar>& massDensity,
                               const FieldSpace::Field<Dimension, Scalar>& temperature) const;

  virtual void setSoundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
                             const FieldSpace::Field<Dimension, Scalar>& massDensity,
                             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setGammaField(FieldSpace::Field<Dimension, Scalar>& gamma,
			     const FieldSpace::Field<Dimension, Scalar>& massDensity,
			     const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setBulkModulus(FieldSpace::Field<Dimension, Scalar>& bulkModulus,
                              const FieldSpace::Field<Dimension, Scalar>& massDensity,
                              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  // Access the member data.
  double a1() const;
  double a2pos() const;
  double a2neg() const;
  double b0() const;
  double b1() const;
  double b2pos() const;
  double b2neg() const;
  double c0() const;
  double c1() const;
  double c2pos() const;
  double c2neg() const;
  double E0() const;
  double atomicWeight() const;
  double Cv() const;

  void a1(const double val);
  void a2pos(const double val);
  void a2neg(const double val);
  void b0(const double val);
  void b1(const double val);
  void b2pos(const double val);
  void b2neg(const double val);
  void c0(const double val);
  void c1(const double val);
  void c2pos(const double val);
  void c2neg(const double val);
  void E0(const double val);
  void atomicWeight(const double val);

  // If requested, the user can specify an external pressure to be applied
  // in the pressure calculation.
  double externalPressure() const;
  void externalPressure(const double P);

  // Compute an individual value for DPDrho.
  double DPDrho(const double massDensity,
                const double specificThermalEnergy) const;

  // Equations of state should have a valid test.
  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mA1, mA2pos, mA2neg, mB0, mB1, mB2pos, mB2neg, mC0, mC1, mC2pos, mC2neg, mE0;
  double mAtomicWeight;
  double mCv;
  double mExternalPressure;

  // No default constructor, copying, or assignment.
  OsborneEquationOfState();
  OsborneEquationOfState(const OsborneEquationOfState&);
  OsborneEquationOfState& operator=(const OsborneEquationOfState&);

  using Material::EquationOfState<Dimension>::mConstants;
};

}
}

#include "OsborneEquationOfStateInline.hh"

#else
// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class OsborneEquationOfState;
  }
}

#endif
