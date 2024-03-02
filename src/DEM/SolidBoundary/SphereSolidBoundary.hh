//---------------------------------Spheral++----------------------------------//
// SphereSolidBoundary -- cylinder with finite length solid boundary for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_SphereSolidBoundary_hh__
#define __Spheral_SphereSolidBoundary_hh__

#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class SphereSolidBoundary : public SolidBoundaryBase<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;

public:
  //--------------------------- Public Interface ---------------------------//

  SphereSolidBoundary(const Vector& center, 
                      const Scalar  radius,
                      const Vector& clipPoint,
                      const Vector& clipAxis);

  ~SphereSolidBoundary();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector velocity(const Vector& position) const override;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;


  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& center() const;
  void center(const Vector& value);

  Scalar radius() const;
  void radius(Scalar value);

  const Vector& clipPoint() const;
  void clipPoint(const Vector& value);

  const Vector& clipAxis() const;
  void clipAxis(const Vector& value);

  const Vector& velocity() const;
  void velocity(const Vector& value);

  void setClipIntersectionRadius();
protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mCenter;
  Scalar mRadius;
  Vector mClipPoint;
  Vector mClipAxis;
  Scalar mClipIntersectionRadius;
  
  Vector mVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  SphereSolidBoundary();
  SphereSolidBoundary(const SphereSolidBoundary&);
  SphereSolidBoundary& operator=(const SphereSolidBoundary&);
};

}

#include "SphereSolidBoundaryInline.hh"

#endif
