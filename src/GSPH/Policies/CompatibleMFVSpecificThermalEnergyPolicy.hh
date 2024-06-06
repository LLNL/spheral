//---------------------------------Spheral++----------------------------------//
// CompatibleMFVSpecificThermalEnergyPolicy -- This is a generalization of the 
//     Lagrangian compatible energy scheme to ALE-based scheme with mass flux
//     between nodes. 
//
// J.M. Pearl 2024
//----------------------------------------------------------------------------//

#ifndef __Spheral_CompatibleMFVSpecificThermalEnergyPolicy_hh__
#define __Spheral_CompatibleMFVSpecificThermalEnergyPolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class CompatibleMFVSpecificThermalEnergyPolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using KeyType = typename UpdatePolicyBase<Dimension>::KeyType;

  // Constructors, destructor.
  CompatibleMFVSpecificThermalEnergyPolicy(const DataBase<Dimension>& db);
  virtual ~CompatibleMFVSpecificThermalEnergyPolicy();
  
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

  CompatibleMFVSpecificThermalEnergyPolicy(const CompatibleMFVSpecificThermalEnergyPolicy& rhs);
  CompatibleMFVSpecificThermalEnergyPolicy& operator=(const CompatibleMFVSpecificThermalEnergyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class CompatibleMFVSpecificThermalEnergyPolicy;
}

#endif
