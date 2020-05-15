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

#include "DataBase/IncrementFieldList.hh"
#include "Geometry/Dimension.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;

class RZNonSymmetricSpecificThermalEnergyPolicy: 
    public IncrementFieldList<Dim<2>, Dim<2>::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::SymTensor SymTensor;
  typedef FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  RZNonSymmetricSpecificThermalEnergyPolicy(const DataBase<Dimension>& db);
  virtual ~RZNonSymmetricSpecificThermalEnergyPolicy();
  
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
