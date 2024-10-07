//---------------------------------Spheral++----------------------------------//
// IncrementASPHHtensor
//
// Specialized version of FieldUpdatePolicy for time integrating the H tensor.
//
// Created by JMO, Mon Oct  7 13:31:02 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementASPHHtensor_hh__
#define __Spheral_IncrementASPHHtensor_hh__

#include "DataBase/FieldUpdatePolicy.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class IncrementASPHHtensor: public FieldUpdatePolicy<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using KeyType = typename FieldUpdatePolicy<Dimension>::KeyType;
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor.
  IncrementASPHHtensor(const Scalar hmin,
                       const Scalar hmax,
                       const Scalar hminratio,
                       const bool fixShape,
                       const bool radialOnly);
  virtual ~IncrementASPHHtensor()   {}
  IncrementASPHHtensor(const IncrementASPHHtensor& rhs) = delete;
  IncrementASPHHtensor& operator=(const IncrementASPHHtensor& rhs) = delete;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt) override;

  // Access the min and max's.
  Scalar hmin() const               { return mhmin; }
  Scalar hmax() const               { return mhmax; }
  Scalar hminratio() const          { return mhminratio; }
  bool fixShape() const             { return mFixShape; }
  bool radialOnly() const           { return mRadialOnly; }

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const override;

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mhmin, mhmax, mhminratio;
  bool mFixShape, mRadialOnly;
};

}

#endif
