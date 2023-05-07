//---------------------------------Spheral++----------------------------------//
// RectangularFinitePlane -- solid planar boundary for DEM with finite extent
//                           and rectangular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_RectangularFinitePlane_hh__
#define __Spheral_RectangularFinitePlane_hh__

#include "DEM/SolidBoundary/SolidBoundary.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class RectangularFinitePlane : public SolidBoundary<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;

public:
  //--------------------------- Public Interface ---------------------------//

  RectangularFinitePlane(const Vector& point,
                         const Vector& exent, 
                         const Tensor& basis);

  ~RectangularFinitePlane();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector velocity(const Vector& position) const override;

  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& point() const;
  void point(const Vector& value);

  const Tensor& basis() const;
  void basis(const Tensor& value);

  const Vector& extent() const;
  void extent(const Vector& value);

  const Vector& velocity() const;
  void velocity(const Vector& value);

protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mPoint;
  Tensor mBasis;
  Vector mExtent;
  Vector mVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  RectangularFinitePlane();
  RectangularFinitePlane(const RectangularFinitePlane&);
  RectangularFinitePlane& operator=(const RectangularFinitePlane&);
};

}

#include "RectangularFinitePlaneInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class RectangularFinitePlane;
}

#endif
