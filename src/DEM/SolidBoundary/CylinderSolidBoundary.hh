//---------------------------------Spheral++----------------------------------//
// CylinderSolidBoundary -- cylinder with finite length solid boundary for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_CylinderSolidBoundary_hh__
#define __Spheral_CylinderSolidBoundary_hh__

#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class CylinderSolidBoundary : public SolidBoundaryBase<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;

public:
  //--------------------------- Public Interface ---------------------------//

  CylinderSolidBoundary(const Vector& point,
                 const Vector& axis, 
                 const Scalar  radius,
                 const Scalar  length);

  ~CylinderSolidBoundary();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector localVelocity(const Vector& position) const override;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

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
  CylinderSolidBoundary();
  CylinderSolidBoundary(const CylinderSolidBoundary&);
  CylinderSolidBoundary& operator=(const CylinderSolidBoundary&);
};

}

#include "CylinderSolidBoundaryInline.hh"

#endif
