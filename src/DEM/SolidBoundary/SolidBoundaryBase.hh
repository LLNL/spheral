//---------------------------------Spheral++----------------------------------//
// SolidBoundaryBase -- A base class for solid wall contact boundaries in DEM.
//                      Solid boundaries are analytically defined surfaces
//                      that behave like soft-sphere contacts. This type of
//                      boundary is distinct from standard spheral boundaries
//                      in that there are no ghost particles. The evolution of 
//                      particle-boundary contact state (i.e sliding/rolling/
//                      torsion springs) is handled exactly like 
//                      particle-particle contacts and the machinery for this
//                      is implemented in the DEMBase class.  The boundaries
//                      themselves simply define the surface and how it evolves
//                      in time.                  
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
