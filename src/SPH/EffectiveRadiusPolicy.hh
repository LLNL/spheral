//---------------------------------Spheral++----------------------------------//
// EffectiveRadiusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the effective radius for RZ coordinates
//
// Created by JMO, Tue Oct 30 15:42:41 PDT 2007
//----------------------------------------------------------------------------//
#ifndef __Spheral_EffectiveRadiusPolicy_hh__
#define __Spheral_EffectiveRadiusPolicy_hh__

#include <string>

#include "DataBase/FieldListUpdatePolicyBase.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;

class EffectiveRadiusPolicy: public FieldListUpdatePolicyBase<Dim<2>, Dim<2>::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;
  using KeyType = FieldListUpdatePolicyBase<Dimension, Scalar>::KeyType;

  // Constructors, destructor.
  EffectiveRadiusPolicy();
  virtual ~EffectiveRadiusPolicy();
  
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
  EffectiveRadiusPolicy(const EffectiveRadiusPolicy& rhs);
  EffectiveRadiusPolicy& operator=(const EffectiveRadiusPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  class EffectiveRadiusPolicy;
}

#endif
