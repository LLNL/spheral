//---------------------------------Spheral++----------------------------------//
// RZNonSymmetricSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy
// as a dependent quantity.
// 
// This version is for use with RZ axisymmetric symmetry.
//
// Created by JMO, Mon Feb  3 21:23:11 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_RZNonSymmetricSpecificThermalEnergyPolicy_hh__
#define __Spheral_RZNonSymmetricSpecificThermalEnergyPolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"
#include "Geometry/Dimension.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;

class RZNonSymmetricSpecificThermalEnergyPolicy: public UpdatePolicyBase<Dim<2>> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using SymTensor = Dimension::SymTensor;
  using KeyType = UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  RZNonSymmetricSpecificThermalEnergyPolicy(const DataBase<Dimension>& db);
  virtual ~RZNonSymmetricSpecificThermalEnergyPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // If the derivative stored values for the pair-accelerations has not been updated,
  // we need to just time advance normally.
  virtual void updateAsIncrement(const KeyType& key,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs,
                                 const double multiplier,
                                 const double t,
                                 const double dt) override;

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  const DataBase<Dimension>* mDataBasePtr;

  RZNonSymmetricSpecificThermalEnergyPolicy(const RZNonSymmetricSpecificThermalEnergyPolicy& rhs);
  RZNonSymmetricSpecificThermalEnergyPolicy& operator=(const RZNonSymmetricSpecificThermalEnergyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class RZNonSymmetricSpecificThermalEnergyPolicy;
}

#endif
