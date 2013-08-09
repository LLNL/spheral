//---------------------------------Spheral++----------------------------------//
// SVPHCorrectionsPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the SVPHCorrections in the state.
//
// Created by JMO, Fri Aug  9 15:24:04 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_SVPHCorrectionsPolicy_hh__
#define __Spheral_SVPHCorrectionsPolicy_hh__

#include <string>

#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/DataBase.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template<typename Dimension>
class SVPHCorrectionsPolicy: public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  SVPHCorrectionsPolicy(const DataBaseSpace::DataBase<Dimension>& dataBase,
                        const KernelSpace::TableKernel<Dimension>& kernel);
  virtual ~SVPHCorrectionsPolicy();
  
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
  const DataBaseSpace::DataBase<Dimension>& mDataBase;
  const KernelSpace::TableKernel<Dimension>& mKernel;

  SVPHCorrectionsPolicy();
  SVPHCorrectionsPolicy(const SVPHCorrectionsPolicy& rhs);
  SVPHCorrectionsPolicy& operator=(const SVPHCorrectionsPolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SVPHCorrectionsPolicy;
}

#endif
