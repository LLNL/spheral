//---------------------------------Spheral++----------------------------------//
// NonSymmetricSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy
// as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method in the case where we do *not* assume that pairwise forces are
// equal and opposite.
//
// We also make this method entirely inlined so the user can specify a special
// weighting function for determining the effective mass.  This is useful for
// curvilinear coordinates.
//
// Created by JMO, Sat Aug 10 23:03:39 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_NonSymmetricSpecificThermalEnergyPolicy_hh__
#define __Spheral_NonSymmetricSpecificThermalEnergyPolicy_hh__

#include "DataBase/IncrementFieldList.hh"

#include <string>
#include <functional>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class NonSymmetricSpecificThermalEnergyPolicy: 
    public IncrementFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  NonSymmetricSpecificThermalEnergyPolicy(const DataBase<Dimension>& db,
                                          std::function<Scalar(const Scalar&, const Vector&)> effectiveMass = [](const Scalar& mi, const Vector& posi) { return mi; });
  virtual ~NonSymmetricSpecificThermalEnergyPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // If the derivative stored values for the pair-accelerations has not been updated,
  // we need to just time advance normally.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) {
    IncrementFieldList<Dimension, Scalar>::update(key,
                                                  state,
                                                  derivs,
                                                  multiplier,
                                                  t,
                                                  dt);
  }

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  const DataBase<Dimension>* mDataBasePtr;
  const std::function<Scalar(const Scalar&, const Vector&)> mEffectiveMassFunc;

  NonSymmetricSpecificThermalEnergyPolicy(const NonSymmetricSpecificThermalEnergyPolicy& rhs);
  NonSymmetricSpecificThermalEnergyPolicy& operator=(const NonSymmetricSpecificThermalEnergyPolicy& rhs);
};

}

#include "Hydro/NonSymmetricSpecificThermalEnergyPolicyInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class NonSymmetricSpecificThermalEnergyPolicy;
}

#endif
