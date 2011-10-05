//---------------------------------Spheral++----------------------------------//
// MeshWall
//
// Class hierarchy to represent the types of boundaries that can be enforced
// on Meshes.  This is meant to mirror (and encapsulate) the wall idea of 
// Voro++.
//
// Created by JMO, Sat Jul  9 23:48:21 PDT 2011
//----------------------------------------------------------------------------//
// Forward declare voro++ stuff.
class wall;
template<typename T> class container_base;
class radius_mono;
typedef container_base<radius_mono> container;
class neighbor_track;
class neighbor_none;
template<typename T> struct voronoicell_base;
typedef voronoicell_base<neighbor_none> voronoicell;

#include "voro++/container.hh"
#include "voro++/wall.hh"
#include "voro++/cell.hh"
#include "boost/shared_ptr.hpp"

#include "MeshWall.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;

//------------------------------------------------------------------------------
// MeshWall
//------------------------------------------------------------------------------
// Constructors.
template<typename Dimension>
MeshWall<Dimension>::
MeshWall():
  mWallPtrs(),
  mXmin(-numeric_limits<double>::max()*Vector::one),
  mXmax( numeric_limits<double>::max()*Vector::one) {
}

template<typename Dimension>
MeshWall<Dimension>::
MeshWall(typename MeshWall<Dimension>::wall_ptr wallPtr):
  mWallPtrs(),
  mXmin(-numeric_limits<double>::max()*Vector::one),
  mXmax( numeric_limits<double>::max()*Vector::one) {
  mWallPtrs.push_back(wallPtr);
}

template<typename Dimension>
MeshWall<Dimension>::
MeshWall(typename MeshWall<Dimension>::wall_ptr wallPtr1, 
         typename MeshWall<Dimension>::wall_ptr wallPtr2):
  mWallPtrs(),
  mXmin(-numeric_limits<double>::max()*Vector::one),
  mXmax( numeric_limits<double>::max()*Vector::one) {
  mWallPtrs.push_back(wallPtr1);
  mWallPtrs.push_back(wallPtr2);
}
  
// Destructor.
template<typename Dimension>
MeshWall<Dimension>::
~MeshWall() {
}

// Return the underlying Voro++ walls.
template<typename Dimension>
vector<typename MeshWall<Dimension>::wall_ptr>&
MeshWall<Dimension>::
wallPtrs() {
  return mWallPtrs;
}

// Add a wall.
template<typename Dimension>
void
MeshWall<Dimension>::
addWall(typename MeshWall<Dimension>::wall_ptr wallPtr) {
  mWallPtrs.push_back(wallPtr);
}

// Check against a new xmin.
template<typename Dimension>
void
MeshWall<Dimension>::
xminApply(const typename Dimension::Vector& xmin) {
  mXmin = elementWiseMax(mXmin, xmin);
}

// Check against a new xmax.
template<typename Dimension>
void
MeshWall<Dimension>::
xmaxApply(const typename Dimension::Vector& xmax) {
  mXmax = elementWiseMin(mXmax, xmax);
}

//------------------------------------------------------------------------------
// Helper method to build min/max Vectors based on coordinate aligned planes.
//------------------------------------------------------------------------------
inline
Dim<2>::Vector
newMinVector(const Dim<2>::Vector& point,
             const Dim<2>::Vector& nhat) {
  typedef Dim<2>::Vector Vector;
  Vector xmin = -numeric_limits<double>::max() * Vector::one;
  if      (fuzzyEqual(nhat.x(),  1.0, 1.0e-5)) { xmin.x(point.x()); }
  else if (fuzzyEqual(nhat.y(),  1.0, 1.0e-5)) { xmin.y(point.y()); }
  return xmin;
}
             
inline
Dim<3>::Vector
newMinVector(const Dim<3>::Vector& point,
             const Dim<3>::Vector& nhat) {
  typedef Dim<3>::Vector Vector;
  Vector xmin = -numeric_limits<double>::max() * Vector::one;
  if      (fuzzyEqual(nhat.x(),  1.0, 1.0e-5)) { xmin.x(point.x()); }
  else if (fuzzyEqual(nhat.y(),  1.0, 1.0e-5)) { xmin.y(point.y()); }
  else if (fuzzyEqual(nhat.z(),  1.0, 1.0e-5)) { xmin.z(point.z()); }
  return xmin;
}
             
inline
Dim<2>::Vector
newMaxVector(const Dim<2>::Vector& point,
             const Dim<2>::Vector& nhat) {
  typedef Dim<2>::Vector Vector;
  Vector xmax = numeric_limits<double>::max() * Vector::one;
  if      (fuzzyEqual(nhat.x(), -1.0, 1.0e-5)) { xmax.x(point.x()); }
  else if (fuzzyEqual(nhat.y(), -1.0, 1.0e-5)) { xmax.y(point.y()); }
  return xmax;
}
             
inline
Dim<3>::Vector
newMaxVector(const Dim<3>::Vector& point,
             const Dim<3>::Vector& nhat) {
  typedef Dim<3>::Vector Vector;
  Vector xmax = numeric_limits<double>::max() * Vector::one;
  if      (fuzzyEqual(nhat.x(), -1.0, 1.0e-5)) { xmax.x(point.x()); }
  else if (fuzzyEqual(nhat.y(), -1.0, 1.0e-5)) { xmax.y(point.y()); }
  else if (fuzzyEqual(nhat.z(), -1.0, 1.0e-5)) { xmax.z(point.z()); }
  return xmax;
}

//------------------------------------------------------------------------------
// PlanarMeshWall
//------------------------------------------------------------------------------
// Constructors.
template<typename Dimension>
PlanarMeshWall<Dimension>::
PlanarMeshWall(const GeomPlane<Dimension>& plane):
  MeshWall<Dimension>(boost::shared_ptr<wall>(new wall_plane(plane.normal().x(),
                                                             plane.normal().y(),
                                                             plane.normal().z(),
                                                             plane.normal().dot(plane.point())))) {
  const Vector xmin = newMinVector(plane.point(), plane.normal());
  const Vector xmax = newMaxVector(plane.point(), plane.normal());
  this->xminApply(xmin);
  this->xmaxApply(xmax);
}

template<typename Dimension>
PlanarMeshWall<Dimension>::
PlanarMeshWall(const GeomPlane<Dimension>& plane1,
               const GeomPlane<Dimension>& plane2):
  MeshWall<Dimension>(boost::shared_ptr<wall>(new wall_plane(plane1.normal().x(),
                                                             plane1.normal().y(),
                                                             plane1.normal().z(),
                                                             plane1.normal().dot(plane1.point()))),
                      boost::shared_ptr<wall>(new wall_plane(plane2.normal().x(),
                                                             plane2.normal().y(),
                                                             plane2.normal().z(),
                                                             plane2.normal().dot(plane2.point())))) {
  const Vector xmin1 = newMinVector(plane1.point(), plane1.normal());
  const Vector xmax1 = newMaxVector(plane1.point(), plane1.normal());
  this->xminApply(xmin1);
  this->xmaxApply(xmax1);
  const Vector xmin2 = newMinVector(plane2.point(), plane2.normal());
  const Vector xmax2 = newMaxVector(plane2.point(), plane2.normal());
  this->xminApply(xmin2);
  this->xmaxApply(xmax2);
}

// Destructor.
template<typename Dimension>
PlanarMeshWall<Dimension>::
~PlanarMeshWall() {
}

//------------------------------------------------------------------------------
// FactedWall
// An implementation of the Voro++ wall for a facted bounding volume.
//------------------------------------------------------------------------------
template<typename FacetedVolume>
class FacetedWall: public wall {
public:
  typedef typename FacetedVolume::Vector Vector;
  typedef typename FacetedVolume::Facet Facet;

  FacetedWall(const FacetedVolume& poly, int id=-99):
    mFacetedVolume(poly),
    mID(id) {}
  virtual ~FacetedWall() {}

  // Test if the given point is inside the wall.
  virtual bool point_inside(fpoint x, fpoint y, fpoint z) {
    return mFacetedVolume.contains(Vector(x, y, z));
  }

  // Cut a cell.
  template<typename Cell>
  bool cut_cell_base(Cell& c, fpoint x, fpoint y, fpoint z) {
    const Vector p(x, y, z);
    const Vector cp = mFacetedVolume.closestPoint(p);
    const Vector delta = cp - p;
    return c.nplane(delta.x(), delta.y(), delta.z(), delta.magnitude2(), mID);
  }

  // The required cut cell methods.
  virtual bool cut_cell(voronoicell_base<neighbor_none>& c, fpoint x, fpoint y, fpoint z) {
    return cut_cell_base(c, x, y, z);
  }

  virtual bool cut_cell(voronoicell_base<neighbor_track>& c, fpoint x, fpoint y, fpoint z) {
    return cut_cell_base(c, x, y, z);
  }

private:
  FacetedVolume mFacetedVolume;
  int mID;
};

//------------------------------------------------------------------------------
// FacetedMeshWall
//------------------------------------------------------------------------------
// Constructors.
template<typename Dimension>
FacetedMeshWall<Dimension>::
FacetedMeshWall(const typename Dimension::FacetedVolume& volume):
  MeshWall<Dimension>(boost::shared_ptr<wall>(new FacetedWall<FacetedVolume>(volume))) {
  this->xminApply(volume.xmin());
  this->xmaxApply(volume.xmax());
}


// Destructor.
template<typename Dimension>
FacetedMeshWall<Dimension>::
~FacetedMeshWall() {
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class MeshWall<Dim<2> >;
template class MeshWall<Dim<3> >;

template class PlanarMeshWall<Dim<2> >;
template class PlanarMeshWall<Dim<3> >;

template class FacetedMeshWall<Dim<2> >;
template class FacetedMeshWall<Dim<3> >;

}
}
