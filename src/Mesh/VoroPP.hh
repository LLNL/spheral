//---------------------------------Spheral++----------------------------------//
// VoroPP
// 
// An interface/container class to encapsulate our use of the Voro++ library.
// This is necessary 'cause Voro++ doesn't build itself, you're supposed to 
// include it's .cc file into your own file.  I can't do this multiple times
// everywhere I use Voro++ because that will cause mutiple instantiations, so 
// hide everything behind a single instantiation in this object.
//
// This header also include the PolyhedralWall object, implementing Voro++'s 
// wall object using our GeomPolyhedron class.
//
// Created by JMO, Mon Apr  4 16:45:17 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_VoroPP__
#define __Spheral_VoroPP__

#include <vector>
#include "boost/unordered_map.hpp"
#include "boost/shared_ptr.hpp"

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "MeshWall.hh"
#include "Cell.hh"

// Forward declare voro++ stuff.
template<typename T> class container_base;
class radius_mono;
typedef container_base<radius_mono> container;
class neighbor_none;
template<typename T> struct voronoicell_base;
typedef voronoicell_base<neighbor_none> voronoicell;
class wall;

namespace Spheral {
namespace MeshSpace {

//------------------------------------------------------------------------------
// VoroPP -- encapsulate the Voro++ container interface.
//------------------------------------------------------------------------------
template<typename Dimension>
class VoroPP {
  //--------------------------- Public Interface ---------------------------//
public:
  typedef typename Dimension::Vector Vector;
  typedef Dim<3>::Vector Vector3d;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  typedef GeomPlane<Dimension> Plane;

  // Constructors, destructor.
  VoroPP(const std::vector<Vector>& generators, const Vector& xmin, const Vector& xmax,
         const unsigned nx = 20,
         const unsigned ny = 20,
         const unsigned nz = 20,
         const double edgeTol = 1.0e-8);
  ~VoroPP();

  // The internal scaling.
  double scale() const;

  // Add a wall.
  void addBoundary(MeshWall<Dimension>& meshWall);

  // Return the vertex and face stuff for all cells.
  // We use our local Cell class to encapsulate the cell by cell info.
  std::vector<unsigned>
  allCells(std::vector<Cell<Dimension> >& cells) const;

  // Compute which sub-region the given position is in.
  void subRegion(const Vector& p, unsigned& i, unsigned& j, unsigned& k) const;

private:
  unsigned mNumGenerators, mNx, mNy, mNz;
  double mEdgeTol, mScale;
  const std::vector<Vector>* mGeneratorsPtr;
  Vector mXmin, mXmax;
  boost::shared_ptr<container> mContainerPtr;

  // Internal method for common constructor operations.
  void construct(const std::vector<Vector>& generators,
                 const unsigned nx,
                 const unsigned ny,
                 const unsigned nz);

  // Compute the appropriate scale factor.
  double computeScale(const Vector& xmin,
                      const Vector& xmax) const;

  // Forbidden methods.
  VoroPP();
  VoroPP(const VoroPP&);
  VoroPP& operator=(const VoroPP&);
};

// Declare the specializations.
template<> VoroPP<Dim<2> >::VoroPP(const std::vector<Dim<2>::Vector>& generators,
                                   const Dim<2>::Vector& xmin,
                                   const Dim<2>::Vector& xmax,
                                   const unsigned nx,
                                   const unsigned ny,
                                   const unsigned nz,
                                   const double edgeTol);
template<> VoroPP<Dim<3> >::VoroPP(const std::vector<Dim<3>::Vector>& generators,
                                   const Dim<3>::Vector& xmin,
                                   const Dim<3>::Vector& xmax,
                                   const unsigned nx,
                                   const unsigned ny,
                                   const unsigned nz,
                                   const double edgeTol);
template<> std::vector<unsigned> VoroPP<Dim<2> >::allCells(std::vector<Cell<Dim<2> > >& cells) const;
template<> std::vector<unsigned> VoroPP<Dim<3> >::allCells(std::vector<Cell<Dim<3> > >& cells) const;

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace MeshSpace {
    template<typename Dimension> class VoroPP;
  }
}

#endif
