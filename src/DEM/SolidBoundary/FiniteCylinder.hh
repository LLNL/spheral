//---------------------------------Spheral++----------------------------------//
// FiniteCylinder -- cylinder with finite length solid boundary for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_FiniteCylinder_hh__
#define __Spheral_FiniteCylinder_hh__

#include "DEM/SolidBoundary/SolidBoundary.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class FiniteCylinder : public SolidBoundary<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;

public:
  //--------------------------- Public Interface ---------------------------//

  FiniteCylinder(const Vector& point,
                 const Vector& axis, 
                 const Scalar  radius,
                 const Scalar  length);

  ~FiniteCylinder();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector velocity(const Vector& position) const override;

  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& point() const;
  void point(const Vector& value);

  const Vector& axis() const;
  void axis(const Vector& value);

  Scalar length() const;
  void length(Scalar value);

  Scalar radius() const;
  void radius(Scalar value);

  const Vector& velocity() const;
  void velocity(const Vector& value);

protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mPoint;
  Vector mAxis;
  Scalar mRadius;
  Scalar mLength;

  Vector mVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  FiniteCylinder();
  FiniteCylinder(const FiniteCylinder&);
  FiniteCylinder& operator=(const FiniteCylinder&);
};

}

#include "FiniteCylinderInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class FiniteCylinder;
}

#endif
