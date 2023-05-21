//---------------------------------Spheral++----------------------------------//
// InfinitePlaneSolidBoundary -- rigid planar wall contact for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_InfinitePlaneSolidBoundary_hh__
#define __Spheral_InfinitePlaneSolidBoundary_hh__

#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

namespace Spheral {


template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class InfinitePlaneSolidBoundary : public SolidBoundaryBase<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;

public:
  //--------------------------- Public Interface ---------------------------//

  InfinitePlaneSolidBoundary(const Vector& point, 
                const Vector& normal);

  ~InfinitePlaneSolidBoundary();

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
  InfinitePlaneSolidBoundary();
  InfinitePlaneSolidBoundary(const InfinitePlaneSolidBoundary&);
  InfinitePlaneSolidBoundary& operator=(const InfinitePlaneSolidBoundary&);
};

}

#include "InfinitePlaneSolidBoundaryInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class InfinitePlaneSolidBoundary;
}

#endif
