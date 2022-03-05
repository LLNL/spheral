//---------------------------------Spheral++----------------------------------//
// CompatibleDifferenceSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy 
// as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method.
//----------------------------------------------------------------------------//

#ifndef __Spheral_CompatibleDifferenceSpecificThermalEnergyPolicy_hh__
#define __Spheral_CompatibleDifferenceSpecificThermalEnergyPolicy_hh__

#include "DataBase/IncrementFieldList.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class CompatibleDifferenceSpecificThermalEnergyPolicy: 
    public IncrementFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  CompatibleDifferenceSpecificThermalEnergyPolicy(const DataBase<Dimension>& db);
  virtual ~CompatibleDifferenceSpecificThermalEnergyPolicy();
  
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

  CompatibleDifferenceSpecificThermalEnergyPolicy(const CompatibleDifferenceSpecificThermalEnergyPolicy& rhs);
  CompatibleDifferenceSpecificThermalEnergyPolicy& operator=(const CompatibleDifferenceSpecificThermalEnergyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class CompatibleDifferenceSpecificThermalEnergyPolicy;
}

#endif
