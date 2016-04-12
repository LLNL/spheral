//---------------------------------Spheral++----------------------------------//
// ANEOS -- An interface to the ANEOS equation of state from 
// Melosh amony others.  This is a C++ wrapper around calls to the underlying
// ANEOS fortran library.  The user must provide the ancillary file containing
// parameters for ANEOS in its expected format.  Also make sure to keep units
// consistent in this separate input file with what you give the C++ EOS here!
//
// Created by JMO, Tue Apr 23 14:55:28 PDT 2013
//----------------------------------------------------------------------------//
#ifndef ANEOS_HH
#define ANEOS_HH

#include "boost/multi_array.hpp"

#include "Material/EquationOfState.hh"

// Forward declarations.
namespace Spheral {
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }
}

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class ANEOS: public Material::EquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename boost::multi_array<double, 2> array_type;
  typedef typename array_type::array_view<1>::type slice_type;
  typedef typename array_type::const_array_view<1>::type const_slice_type;
  typedef boost::multi_array_types::index_range range;

  // Constructors, destructors.
  ANEOS(const int materialNumber,
        const unsigned numRhoVals,
        const unsigned numTvals,
        const double rhoMin,
        const double rhoMax,
        const double Tmin,
        const double Tmax,
        const Material::PhysicalConstants& constants,
        const double externalPressure,
        const double minimumPressure,
        const double maximumPressure,
        const Material::MaterialPressureMinType minPressureType);
  ~ANEOS();

  // We require any equation of state to define the following properties.
  virtual void setPressure(FieldSpace::Field<Dimension, Scalar>& Pressure,
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

  virtual void setEntropy(FieldSpace::Field<Dimension, Scalar>& entropy,
                          const FieldSpace::Field<Dimension, Scalar>& massDensity,
                          const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  // We also want the equivalent functions for individual calculations.
  Scalar pressure(const Scalar massDensity,
                  const Scalar specificThermalEnergy) const;

  Scalar temperature(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar specificThermalEnergy(const Scalar massDensity,
                               const Scalar temperature) const;

  Scalar specificHeat(const Scalar massDensity,
                      const Scalar temperature) const;

  Scalar soundSpeed(const Scalar massDensity,
                    const Scalar specificThermalEnergy) const;

  // Get the effective gamma (ratio of specific heats) for this eos.
  Scalar gamma(const Scalar massDensity,
               const Scalar specificThermalEnergy) const;

  // Get the bulk modulus.
  Scalar bulkModulus(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar entropy(const Scalar massDensity,
                 const Scalar specificThermalEnergy) const;

  // The valid method.
  virtual bool valid() const;

  // Access local variables used to lookup eps based on T.
  int materialNumber() const;
  unsigned numRhoVals() const;
  unsigned numTvals() const;
  double rhoMin() const;
  double rhoMax() const;
  double Tmin() const;
  double Tmax() const;
  const array_type& specificThermalEnergyVals() const;

  // If requested, the user can specify an external pressure to be applied
  // in the pressure calculation.
  double externalPressure() const;
  void externalPressure(const double x);

  double atomicWeight() const;

private:
  //--------------------------- Private Interface ---------------------------//
  // Tables for the temp->energy lookup.
  int mMaterialNumber;
  unsigned mNumRhoVals, mNumTvals;
  double mRhoMin, mRhoMax, mTmin, mTmax, mExternalPressure;
  array_type mSTEvals;

  // ANEOS internal units.
  Material::PhysicalConstants mANEOSunits;

  // Units conversion from ANEOS.
  double mRhoConv, mTconv, mPconv, mEconv, mCVconv, mVelConv, mSconv;

  // Atomic weight.
  double mAtomicWeight;

  // Disallow default constructor
  ANEOS();

  using Material::EquationOfState<Dimension>::mConstants;
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class ANEOS;
  }
}

#endif
