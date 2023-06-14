//---------------------------------Spheral++----------------------------------//
// SolidBoundaryBase -- A base class for solid wall contact bcs
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_SolidBoundaryBase_hh__
#define __Spheral_SolidBoundaryBase_hh__

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class SolidBoundaryBase {

public:
//--------------------------- Public Interface ---------------------------//

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  SolidBoundaryBase();

  virtual ~SolidBoundaryBase();

  virtual Vector distance(const Vector& position) const = 0;
  virtual Vector velocity(const Vector& position) const = 0;

  virtual void update(const double multiplier,
                      const double t,
                      const double dt) = 0;
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SolidBoundaryBase;
}

#endif
