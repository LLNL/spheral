//---------------------------------Spheral++----------------------------------//
// SpecificFromTotalThermalEnergyPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the specific thermal energy from the total
// energy.
// 
// Created by JMO, Thu Jan  7 23:08:39 PST 2016
//----------------------------------------------------------------------------//

#ifndef __Spheral_SpecificFromTotalThermalEnergyPolicy_hh__
#define __Spheral_SpecificFromTotalThermalEnergyPolicy_hh__

#include <string>

#include "DataBase/FieldListUpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
namespace NodeSpace {
  template<typename Dimension> class FluidNodeList;
}
namespace FieldSpace {
  template<typename Dimension, typename DataType> class FieldList;
}
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}

template<typename Dimension>
class SpecificFromTotalThermalEnergyPolicy: 
    public FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  SpecificFromTotalThermalEnergyPolicy();
  virtual ~SpecificFromTotalThermalEnergyPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  SpecificFromTotalThermalEnergyPolicy(const SpecificFromTotalThermalEnergyPolicy& rhs);
  SpecificFromTotalThermalEnergyPolicy& operator=(const SpecificFromTotalThermalEnergyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SpecificFromTotalThermalEnergyPolicy;
}

#endif
