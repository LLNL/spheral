//---------------------------------Spheral++----------------------------------//
// IncrementCullenMultipliers
//
// Created by JMO, Mon Dec 28 22:18:44 PST 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementCullenMultipliers_hh__
#define __Spheral_IncrementCullenMultipliers_hh__

#include <float.h>
#include "DataBase/IncrementBoundedFieldList.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class IncrementCullenMultipliers: public IncrementBoundedFieldList<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  IncrementCullenMultipliers(const Scalar minValue, const Scalar maxValue, const bool hopkinsForm);
  virtual ~IncrementCullenMultipliers();
  
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
  bool mHopkinsForm;
  IncrementCullenMultipliers(const IncrementCullenMultipliers& rhs);
  IncrementCullenMultipliers& operator=(const IncrementCullenMultipliers& rhs);
};

}

#endif
