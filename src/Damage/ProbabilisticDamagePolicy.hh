//---------------------------------Spheral++----------------------------------//
// ProbabilisticDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent scalar damage state.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//   Randles, P.W. & Libersky, L.D., 1996, "Comput. Methods Appl. Engrg, 
//     139, 375-408
//
// Created by JMO, Sun Sep 26 16:15:17 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ProbabilisticDamagePolicy_hh__
#define __Spheral_ProbabilisticDamagePolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class ProbabilisticDamagePolicy: 
    public UpdatePolicyBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename UpdatePolicyBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  explicit ProbabilisticDamagePolicy(const bool damageInCompression,  // allow damage in compression
                                     const double kWeibull,           // coefficient in Weibull power-law
                                     const double mWeibull,           // exponenent in Weibull power-law
                                     const size_t minFlawsPerNode,    // minimum number of flaws to seed on any node
                                     const double Vmin,               // minimum (initial) node volume
                                     const double Vmax);              // maximum (initial) node volume
  virtual ~ProbabilisticDamagePolicy();
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return "delta "; }

private:
  //--------------------------- Private Interface ---------------------------//
  bool mDamageInCompression;
  size_t mMinFlawsPerNode;
  double mkWeibull, mmWeibull, mVmin, mVmax;

  ProbabilisticDamagePolicy();
  ProbabilisticDamagePolicy(const ProbabilisticDamagePolicy& rhs);
  ProbabilisticDamagePolicy& operator=(const ProbabilisticDamagePolicy& rhs);
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class ProbabilisticDamagePolicy;
}

#endif
