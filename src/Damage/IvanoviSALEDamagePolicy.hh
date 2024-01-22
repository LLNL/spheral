//---------------------------------Spheral++----------------------------------//
// IvanoviSALEDamageModel
//
// The Ivanov damage model, hopefully close to how it's implemented in iSALE.
// This damage model is most appropriate for rocky materials.
//
// Refs:
// 
// Collins, G. S., Melosh, H. J., & Ivanov, B. A. (2004). Modeling damage and deformation in impact simulations.
//   Meteoritics & Planetary Science, 39(2), 217–231. http://doi.wiley.com/10.1111/j.1945-5100.2004.tb00337.x
//
// Raducan, S. D., Davison, T. M., Luther, R., & Collins, G. S. (2019). The role of asteroid strength, porosity and
//   internal friction in impact momentum transfer. Icarus. https://doi.org/10.1016/J.ICARUS.2019.03.040
//
// Lundborg, N. (1967). The strength-size relation of granite. International Journal of Rock Mechanics and
//   Mining Sciences & Geomechanics Abstracts, 4(3):269–272.f
//
//---------------------------------Spheral++----------------------------------//
// IvanoviSALEDamageModel
//
// The Ivanov damage model, hopefully close to how it's implemented in iSALE.
// This damage model is most appropriate for rocky materials.
//
// Refs:
// 
// Collins, G. S., Melosh, H. J., & Ivanov, B. A. (2004). Modeling damage and deformation in impact simulations.
//   Meteoritics & Planetary Science, 39(2), 217–231. http://doi.wiley.com/10.1111/j.1945-5100.2004.tb00337.x
//
// Raducan, S. D., Davison, T. M., Luther, R., & Collins, G. S. (2019). The role of asteroid strength, porosity and
//   internal friction in impact momentum transfer. Icarus. https://doi.org/10.1016/J.ICARUS.2019.03.040
//
// Lundborg, N. (1967). The strength-size relation of granite. International Journal of Rock Mechanics and
//   Mining Sciences & Geomechanics Abstracts, 4(3):269–272.f
//
//----------------------------------------------------------------------------//
#ifndef __Spheral_IvanoviSALEDamagePolicy_hh__
#define __Spheral_IvanoviSALEDamagePolicy_hh__

#include "DataBase/UpdatePolicyBase.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class IvanoviSALEDamagePolicy: 
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
  explicit IvanoviSALEDamagePolicy(const double minPlasticFailure,             // minimum plastic strain for failure
                                   const double plasticFailurePressureSlope,   // slope for critical plastic strain
                                   const double plasticFailurePressureOffset,  // intercept for critical plastic strain
                                   const double tensileFailureStress);         // threshold for tensile failure
  virtual ~IvanoviSALEDamagePolicy();
  
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
  double mEpsPfb, mB, mPc, mTensileFailureStress;

  IvanoviSALEDamagePolicy();
  IvanoviSALEDamagePolicy(const IvanoviSALEDamagePolicy& rhs);
  IvanoviSALEDamagePolicy& operator=(const IvanoviSALEDamagePolicy& rhs);
};

}

#endif
