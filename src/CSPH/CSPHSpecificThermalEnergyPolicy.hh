//---------------------------------Spheral++----------------------------------//
// CSPHSpecificThermalEnergyPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the specific thermal energy as a dependent 
// quantity.
// 
// This version is specialized for the compatible energy discretization 
// for use with CSPH.
//
// Created by JMO, Sat Aug  2 06:44:50 PDT 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_CSPHSpecificThermalEnergyPolicy_hh__
#define __Spheral_CSPHSpecificThermalEnergyPolicy_hh__

#include <string>

#include "DataBase/IncrementFieldList.hh"
#include "Kernel/TableKernel.hh"

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
class CSPHSpecificThermalEnergyPolicy: 
    public IncrementFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  CSPHSpecificThermalEnergyPolicy(const DataBaseSpace::DataBase<Dimension>& db,
                                  const KernelSpace::TableKernel<Dimension>& W);
  virtual ~CSPHSpecificThermalEnergyPolicy();
  
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
  const KernelSpace::TableKernel<Dimension>& mWT;

  CSPHSpecificThermalEnergyPolicy(const CSPHSpecificThermalEnergyPolicy& rhs);
  CSPHSpecificThermalEnergyPolicy& operator=(const CSPHSpecificThermalEnergyPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class CSPHSpecificThermalEnergyPolicy;
}

#endif
