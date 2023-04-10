//---------------------------------Spheral++----------------------------------//
// SolidBoundary -- this is a base class for solid wall contact bcs
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidBoundary_hh__
#define __Spheral_SolidBoundary_hh__

//#include <string>

#include <vector>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class SolidBoundary {

public:
//--------------------------- Public Interface ---------------------------//

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  SolidBoundary();

  virtual ~SolidBoundary();

  virtual Scalar value(const Vector& position) const = 0;

  virtual void update(const double multiplier,
                      const double t,
                      const double dt) = 0;


private:

//  std::vector<SolidBoundary<Dimension>*> mClipSurfaces;
};



}

#include "SolidBoundaryInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SolidBoundary;
}

#endif
