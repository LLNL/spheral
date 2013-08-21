//---------------------------------Spheral++----------------------------------//
// PlanarBoundary -- Abstract base class for the boundary conditions defined
// on planes.  Planar boundaries in general are defined in terms of parallel
// pairs of planes, representing entrance and exit conditions.
//
// Created by JMO, Thu Mar  2 21:34:25 PST 2000
//----------------------------------------------------------------------------//

#ifndef PlanarBoundary_HH
#define PlanarBoundary_HH

#include <string>

#ifndef __GCCXML__
#include <vector>
#include "DataOutput/registerWithRestart.hh"
#else
#include "fakestl.hh"
#endif

#include "Boundary.hh"

namespace Spheral {
  template<typename Dimension> class GeomPlane;
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class PlanarBoundary: public Boundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Boundary<Dimension>::BoundaryNodes BoundaryNodes;

  // Constructors and destructors.
  PlanarBoundary();
  PlanarBoundary(const GeomPlane<Dimension>& enterPlane,
                 const GeomPlane<Dimension>& exitPlane);
  virtual ~PlanarBoundary();

  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);
  virtual void updateGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Find the set of nodes in violation of this boundary in the given NodeList.
  // For planar boundaries this is any node that is "behind" the enter plane.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Set the ghost nodes for a predefined set of control nodes.
  void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList, 
                     const std::vector<int>& presetControlNodes);

  // Allow access to the entrance and exit planes.
  virtual const GeomPlane<Dimension>& enterPlane() const;
  virtual void setEnterPlane(const GeomPlane<Dimension>& enterPlane);

  virtual const GeomPlane<Dimension>& exitPlane() const;
  virtual void setExitPlane(const GeomPlane<Dimension>& exitPlane);

  // Function to map positions through the two planes (from the enter 
  // through the exit plane).
  Vector mapPosition(const Vector& position,
		     const GeomPlane<Dimension>& enterPlane,
		     const GeomPlane<Dimension>& exitPlane) const;

  // Determine if the boundary is in a "valid", ready to use state.
  virtual bool valid() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "PlanarBoundary"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Override the clip method for clipping a box.
  virtual void clip(Vector& xmin, Vector& xmax) const;

  // Provide a method to identify tessellation faces on a plane.
  std::vector<unsigned> facesOnPlane(const MeshSpace::Mesh<Dimension>& mesh,
                                     const GeomPlane<Dimension>& plane,
                                     const Scalar tol) const;

private:
  //--------------------------- Private Interface ---------------------------//
  GeomPlane<Dimension> mEnterPlane;
  GeomPlane<Dimension> mExitPlane;

#ifndef __GCCXML__
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif

  // Method to set the ghost node indicies for a given NodeList once the
  // master nodes are set.
  void setGhostNodeIndicies(NodeSpace::NodeList<Dimension>& nodeList);
};

}
}

#ifndef __GCCXML__
#include "PlanarBoundaryInline.hh"
#endif

#else

namespace Spheral {
  namespace BoundarySpace {
    // Forward declaration.
    template<typename Dimension> class PlanarBoundary;
  }
}

#endif
