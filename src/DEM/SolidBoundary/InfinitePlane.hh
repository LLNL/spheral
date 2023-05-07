//---------------------------------Spheral++----------------------------------//
// InfinitePlane -- rigid planar wall contact for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_InfinitePlane_hh__
#define __Spheral_InfinitePlane_hh__

#include "DEM/SolidBoundary/SolidBoundary.hh"

namespace Spheral {


template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class InfinitePlane : public SolidBoundary<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;

public:
  //--------------------------- Public Interface ---------------------------//

  InfinitePlane(const Vector& point, 
                const Vector& normal);

  ~InfinitePlane();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector velocity(const Vector& position) const override;

  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& point() const;
  void point(const Vector& value);

  const Vector& normal() const;
  void normal(const Vector& value);

  const Vector& velocity() const;
  void velocity(const Vector& value);

protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mPoint;
  Vector mNormal;
  Vector mVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  InfinitePlane();
  InfinitePlane(const InfinitePlane&);
  InfinitePlane& operator=(const InfinitePlane&);
};

}

#include "InfinitePlaneInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class InfinitePlane;
}

#endif
