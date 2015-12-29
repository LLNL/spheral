//---------------------------------Spheral++----------------------------------//
// GenericBodyForce -- The base class for all Spheral++ body force
// implementations.
//
// Created by JMO, Wed May 24 14:23:10 PDT 2000
//----------------------------------------------------------------------------//
#ifndef GenericBodyForce_HH
#define GenericBodyForce_HH

#include "Physics.hh"

namespace Spheral {
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

namespace PhysicsSpace {

template<typename Dimension>
class GenericBodyForce: public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  GenericBodyForce();

  // Destructor.
  virtual ~GenericBodyForce();

  // Provide default methods for creating and registering an acceleration source.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state);
  virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);

  // Access to the derivative fields.
  const FieldSpace::FieldList<Dimension, Vector>& DxDt() const;
  const FieldSpace::FieldList<Dimension, Vector>& DvDt() const;

private:
  //--------------------------- Public Interface ---------------------------//
  // Our derivative fields.
  FieldSpace::FieldList<Dimension, Vector> mDxDt;
  FieldSpace::FieldList<Dimension, Vector> mDvDt;
};
}
}

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class GenericBodyForce;
  }
}

#endif
