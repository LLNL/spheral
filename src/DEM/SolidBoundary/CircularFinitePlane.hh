//---------------------------------Spheral++----------------------------------//
// CircularFinitePlane -- solid planar boundary for DEM with finite extent
//                           and circular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_CircularFinitePlane_hh__
#define __Spheral_CircularFinitePlane_hh__

#include "DEM/SolidBoundary/SolidBoundary.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class CircularFinitePlane : public SolidBoundary<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;

public:
  //--------------------------- Public Interface ---------------------------//

  CircularFinitePlane(const Vector& point,
                      const Vector& normal,
                      const Scalar& exent);

  ~CircularFinitePlane();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector velocity(const Vector& position) const override;

  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& point() const;
  void point(const Vector& value);

    const Vector& normal() const;
  void normal(const Vector& value);

  const Scalar& extent() const;
  void extent(const Scalar& value);

  const Vector& velocity() const;
  void velocity(const Vector& value);

protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mPoint;
  Vector mNormal;
  Scalar mExtent;
  Vector mVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  CircularFinitePlane();
  CircularFinitePlane(const CircularFinitePlane&);
  CircularFinitePlane& operator=(const CircularFinitePlane&);
};

}

#include "CircularFinitePlaneInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class CircularFinitePlane;
}

#endif
