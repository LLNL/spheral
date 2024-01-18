//---------------------------------Spheral++----------------------------------//
// OmegaGradhPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the grad h correction terms.
//
// Created by JMO, Tue Oct 30 15:42:41 PDT 2007
//----------------------------------------------------------------------------//
#ifndef __Spheral_OmegaGradhPolicy_hh__
#define __Spheral_OmegaGradhPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;

template<typename Dimension>
class OmegaGradhPolicy: public UpdatePolicyBase<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Field<Dimension, Scalar> FieldType;
  typedef typename UpdatePolicyBase<Dimension, FieldType>::KeyType KeyType;

  // Constructors, destructor.
  OmegaGradhPolicy(const DataBase<Dimension>& dataBase);
  virtual ~OmegaGradhPolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension, Scalar>& rhs) const;

private:
  //--------------------------- Private Interface ---------------------------//
  const DataBase<Dimension>& mDataBase;

  OmegaGradhPolicy(const OmegaGradhPolicy& rhs);
  OmegaGradhPolicy& operator=(const OmegaGradhPolicy& rhs);
};

}

#endif
