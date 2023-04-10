//---------------------------------Spheral++----------------------------------//
// PlanarWall -- this is a base class for solid wall contact bcs
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_PlanarWall_hh__
#define __Spheral_PlanarWall_hh__

//#include <string>

#include "DEM/SolidBoundary/SolidBoundary.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class PlanarWall : public SolidBoundary<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;

public:
  //--------------------------- Public Interface ---------------------------//

  PlanarWall(const Vector& point, 
             const Vector& normal);

  ~PlanarWall();

  virtual Scalar value(const Vector& position) const override;

  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& point() const;
  void point(const Vector& value);

  const Vector& normal() const;
  void normal(const Vector& value);

protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mPoint;
  Vector mNormal;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  PlanarWall();
  PlanarWall(const PlanarWall&);
  PlanarWall& operator=(const PlanarWall&);
};

}

#include "PlanarWallInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PlanarWall;
}

#endif
