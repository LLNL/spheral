//---------------------------------Spheral++----------------------------------//
// MeshWall
//
// Class hierarchy to represent the types of boundaries that can be enforced
// on Meshes.  This is meant to mirror (and encapsulate) the wall idea of 
// Voro++.
//
// Because 1D, 2D, and 3D are so different in how we generate the meshes,
// for now we're creating individual classes with unique interfaces for each
// of these cases.
//
// Created by JMO, Sat Jul  9 23:48:21 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_MeshWall__
#define __Spheral_MeshWall__

// Forward declare voro++ stuff.
class wall;
template<typename T> class container_base;
class radius_mono;
typedef container_base<radius_mono> container;
class neighbor_track;
class neighbor_none;
template<typename T> struct voronoicell_base;
typedef voronoicell_base<neighbor_none> voronoicell;

#include <vector>
#include <limits>
#include "boost/shared_ptr.hpp"
#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

namespace Spheral {
namespace MeshSpace {

//------------------------------------------------------------------------------
// MeshWall
//------------------------------------------------------------------------------
// 2D and 3D.
template<typename Dimension>
class MeshWall  {
public:
  typedef boost::shared_ptr<wall> wall_ptr;
  typedef typename Dimension::Vector Vector;

  // Constructors, destructors, assignment.
  MeshWall();
  virtual ~MeshWall();

  // In 2D and 3D our only real job is to contain Voro++ wall objects.
  std::vector<wall_ptr>& wallPtrs();
  void addWall(wall_ptr wallPtr);

  // Return the current bounding box.
  const Vector& xmin() const { return mXmin; }
  const Vector& xmax() const { return mXmax; }

protected:
  // Convenient constructors for typical descendents.
  MeshWall(wall_ptr wallPtr);
  MeshWall(wall_ptr wallPtr1,
           wall_ptr wallPtr2);

  void xminApply(const Vector& x);
  void xmaxApply(const Vector& x);

private:
  std::vector<boost::shared_ptr<wall> > mWallPtrs;
  Vector mXmin, mXmax;

  // Forbidden methods.
  MeshWall(const MeshWall& rhs);
  MeshWall& operator=(const MeshWall& rhs);
};

//------------------------------------------------------------------------------
// 1D
template<>
class MeshWall<Dim<1> >  {
public:
  typedef Dim<1>::Vector Vector;

  // Constructors, destructors, assignment.
  MeshWall(const double xmin = -std::numeric_limits<double>::max(),
           const double xmax =  std::numeric_limits<double>::max()):
    mXmin(xmin),
    mXmax(xmax) {};
  virtual ~MeshWall() {};

  // Return the min and max allowed x values.
  const Vector& xmin() { return mXmin; }
  const Vector& xmax() { return mXmax; }

private:
  Vector mXmin, mXmax;

  // Forbidden methods.
  MeshWall(const MeshWall& rhs);
  MeshWall& operator=(const MeshWall& rhs);
};

//------------------------------------------------------------------------------
// PlanarMeshWall
//------------------------------------------------------------------------------
// 2D and 3D.
template<typename Dimension>
class PlanarMeshWall: public MeshWall<Dimension>  {
public:
  typedef GeomPlane<Dimension> Plane;
  typedef typename Dimension::Vector Vector;

  // Constructors, destructors, assignment.
  PlanarMeshWall(const Plane& plane);
  PlanarMeshWall(const Plane& plane1,
                 const Plane& plane2);
  virtual ~PlanarMeshWall();

private:
  // Forbidden methods.
  PlanarMeshWall();
  PlanarMeshWall(const PlanarMeshWall& rhs);
  PlanarMeshWall& operator=(const PlanarMeshWall& rhs);
};

//------------------------------------------------------------------------------
// 1D
template<>
class PlanarMeshWall<Dim<1> >: public MeshWall<Dim<1> >  {
public:
  typedef GeomPlane<Dim<1> > Plane;

  // Constructors, destructors, assignment.
  PlanarMeshWall(const Plane& plane):
    MeshWall<Dim<1> >((plane.normal().x() > 0.0 ? plane.point().x() : -std::numeric_limits<double>::max()),
                      (plane.normal().x() > 0.0 ? std::numeric_limits<double>::max() : plane.point().x())) {};
  PlanarMeshWall(const Plane& plane1,
                 const Plane& plane2):
    MeshWall<Dim<1> >(std::max((plane1.normal().x() > 0.0 ? plane1.point().x() : -std::numeric_limits<double>::max()),
                               (plane2.normal().x() > 0.0 ? plane2.point().x() : -std::numeric_limits<double>::max())),
                      std::min((plane1.normal().x() > 0.0 ? std::numeric_limits<double>::max() : plane1.point().x()),
                               (plane2.normal().x() > 0.0 ? std::numeric_limits<double>::max() : plane2.point().x()))) {}
  virtual ~PlanarMeshWall() {};

private:
  // Forbidden methods.
  PlanarMeshWall();
  PlanarMeshWall(const PlanarMeshWall& rhs);
  PlanarMeshWall& operator=(const PlanarMeshWall& rhs);
};

//------------------------------------------------------------------------------
// FacetedMeshWall
//------------------------------------------------------------------------------
// 2D and 3D.
template<typename Dimension>
class FacetedMeshWall: public MeshWall<Dimension>  {
public:
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Constructors, destructors, assignment.
  FacetedMeshWall(const FacetedVolume& volume);
  virtual ~FacetedMeshWall();

private:
  // Forbidden methods.
  FacetedMeshWall();
  FacetedMeshWall(const FacetedMeshWall& rhs);
  FacetedMeshWall& operator=(const FacetedMeshWall& rhs);
};

//------------------------------------------------------------------------------
// 1D
template<>
class FacetedMeshWall<Dim<1> >: public MeshWall<Dim<1> >  {
public:
  typedef Dim<1>::FacetedVolume FacetedVolume;

  // Constructors, destructors, assignment.
  FacetedMeshWall(const FacetedVolume& volume):
    MeshWall<Dim<1> >(volume.xmin().x(), volume.xmax().x()) {};
  virtual ~FacetedMeshWall() {};

private:
  // Forbidden methods.
  FacetedMeshWall();
  FacetedMeshWall(const FacetedMeshWall& rhs);
  FacetedMeshWall& operator=(const FacetedMeshWall& rhs);
};

}
}

#endif
