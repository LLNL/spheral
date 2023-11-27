//---------------------------------Spheral++----------------------------------//
// PositionPolicy -- An implementation of UpdatePolicyBase specialized for the
// updating the position.
//
// This version ignores the XSPH approximation in order to time center the
// velocity for updating the position.  This is intended for use with the 
// compatible energy evolution hydro approximation.
//
// Created by JMO, Mon Jun 19 22:06:07 PDT 2006
//----------------------------------------------------------------------------//
#ifndef __Spheral_PositionPolicy_hh__
#define __Spheral_PositionPolicy_hh__

#include "DataBase/IncrementFieldList.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class PositionPolicy: 
    public IncrementFieldList<Dimension, typename Dimension::Vector> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename FieldListUpdatePolicyBase<Dimension, Vector>::KeyType KeyType;

  // Constructors, destructor.
  PositionPolicy();
  virtual ~PositionPolicy();
  
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
  PositionPolicy(const PositionPolicy& rhs);
  PositionPolicy& operator=(const PositionPolicy& rhs);
};

}

#endif
