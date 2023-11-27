//---------------------------------Spheral++----------------------------------//
// CellPressurePolicy
//
// Created by JMO, Sun Aug 25 14:40:42 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_CellPressurePolicy_hh__
#define __Spheral_CellPressurePolicy_hh__

#include <string>

#include "DataBase/FieldUpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class CellPressurePolicy: public FieldUpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename FieldUpdatePolicyBase<Dimension, Scalar>::KeyType KeyType;

  // Constructors, destructor.
  CellPressurePolicy();
  virtual ~CellPressurePolicy();
  
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
  CellPressurePolicy(const CellPressurePolicy& rhs);
  CellPressurePolicy& operator=(const CellPressurePolicy& rhs);
};

}

#endif
