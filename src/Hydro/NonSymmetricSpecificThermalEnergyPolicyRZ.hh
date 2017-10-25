//---------------------------------Spheral++----------------------------------//
// NonSymmetricSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy
// as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method in RZ coordinates, which implicitly means we have to use the
// non-symmetric (pairwise) form.
//
// Created by JMO, Wed May  4 16:49:59 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_NonSymmetricSpecificThermalEnergyPolicyRZ_hh__
#define __Spheral_NonSymmetricSpecificThermalEnergyPolicyRZ_hh__

#include <string>

#include "Geometry/Dimension.hh"
#include "DataBase/IncrementFieldList.hh"

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

class NonSymmetricSpecificThermalEnergyPolicyRZ: 
    public IncrementFieldList<Dim<2>, Dim<2>::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  NonSymmetricSpecificThermalEnergyPolicyRZ(const DataBaseSpace::DataBase<Dimension>& db);
  virtual ~NonSymmetricSpecificThermalEnergyPolicyRZ();
  
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
  const DataBaseSpace::DataBase<Dimension>* mDataBasePtr;

  NonSymmetricSpecificThermalEnergyPolicyRZ(const NonSymmetricSpecificThermalEnergyPolicyRZ& rhs);
  NonSymmetricSpecificThermalEnergyPolicyRZ& operator=(const NonSymmetricSpecificThermalEnergyPolicyRZ& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  class NonSymmetricSpecificThermalEnergyPolicyRZ;
}

#endif
