//---------------------------------Spheral++----------------------------------//
// CRKSPHVoidBoundary -- Special boundary class to create void points off of
// the surface.
//----------------------------------------------------------------------------//
#include "CRKSPHVoidBoundary.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/FieldBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {
namespace BoundarySpace {

using namespace std;
using std::vector;

using NodeSpace::NodeList;
using FieldSpace::FieldBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHVoidBoundary<Dimension>::
CRKSPHVoidBoundary(const FieldList<Dimension, int>& surfacePoint,
                   const FieldList<Dimension, Vector>& m1):
  Boundary<Dimension>(),
  mSurfacePoint(surfacePoint),
  mM1(m1) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHVoidBoundary<Dimension>::~CRKSPHVoidBoundary() {
}

//------------------------------------------------------------------------------
// setGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
setGhostNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);

  typename Boundary<Dimension>::BoundaryNodes& boundaryNodes = this->accessBoundaryNodes(nodeList);
  vector<int>& cNodes = boundaryNodes.controlNodes;
  vector<int>& gNodes = boundaryNodes.ghostNodes;
  cNodes = vector<int>();
  gNodes = vector<int>();
  const unsigned firstNewGhostNode = nodeList.numNodes();
    
  // Line up the surface points as controls for the void points.
  unsigned numSurfacePoints = 0;
  unsigned j = firstNewGhostNode;
  if (mSurfacePoint.haveNodeList(nodeList)) {
    const Field<Dimension, int>& surfacePoint = **mSurfacePoint.fieldForNodeList(nodeList);
    for (unsigned i = 0; i < firstNewGhostNode; ++i) {
      if (surfacePoint(i) == 1) {
        ++numSurfacePoints;
        cNodes.push_back(i);
        gNodes.push_back(j++);
      }
    }
  }
  CHECK(cNodes.size() == numSurfacePoints);
  CHECK(gNodes.size() == numSurfacePoints);

  // Create the correct number of new void ghost points.
  const unsigned numGhostNodes0 = nodeList.numGhostNodes();
  nodeList.numGhostNodes(numGhostNodes0 + numSurfacePoints);

  // Update the void point positions and H's.
  this->updateGhostNodes(nodeList);
}

//------------------------------------------------------------------------------
// updateGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
updateGhostNodes(NodeList<Dimension>& nodeList) {
  const vector<int>& cNodes = this->controlNodes(nodeList);
  const vector<int>& gNodes = this->ghostNodes(nodeList);
  const unsigned nsurf = cNodes.size();
  CHECK(gNodes.size() == nsurf);
  
  if (nsurf > 0) {
    const Scalar nPerh = nodeList.nodesPerSmoothingScale();
    const Field<Dimension, Vector>& m1 = **mM1.fieldForNodeList(nodeList);
    Field<Dimension, Vector>& pos = nodeList.positions();
    Field<Dimension, SymTensor>& H = nodeList.Hfield();

    // Update the void points positions based on projecting from the surface point.
    for (unsigned k = 0; k != nsurf; ++k) {
      const unsigned i = cNodes[k];
      const Vector& xi = pos(i);
      const SymTensor& Hi = H(i);
      const Vector nhat = m1(i).unitVector();

      const unsigned j = gNodes[k];
      Vector& xj = pos(j);
      SymTensor& Hj = H(j);

      xj = xi + Hi.Inverse()*nhat/nPerh;
      Hj = Hi;
      // cerr << "Void : " << i << " -> " << j << "     " << xi << " -> " << xj << endl;
    }
  }
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, int>& field) const {
  const vector<int>& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  if (field.name() == HydroFieldNames::voidPoint) {
    // voidPoint: flag only ghost void points
    for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = 1;
  } else {
    for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = 0;
  }
}

// Specialization for scalar fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
  const vector<int>& cNodes = this->controlNodes(field.nodeList());
  const vector<int>& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nsurf = cNodes.size();
  CHECK(gNodes.size() == nsurf);
  if (field.name() == HydroFieldNames::volume) {
    // Volume: copy surface->void
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = field(cNodes[k]);
  } else if (field.name() == HydroFieldNames::mass or
             field.name() == HydroFieldNames::massDensity) {
    // mass, mass density: negligible but non-zero
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = std::numeric_limits<Scalar>::epsilon();
  } else {
    // Default zero.
    for (unsigned k = 0; k < nsurf; ++k) field(gNodes[k]) = 0.0;
  }
}

// Specialization for Vector fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  const vector<int>& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = Vector::zero;
}

// Specialization for Tensor fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  const vector<int>& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = Tensor::zero;
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  const vector<int>& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = SymTensor::zero;
}

// Specialization for third rank tensors.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  const vector<int>& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = ThirdRankTensor::zero;
}

// Specialization for vector<scalar> fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, std::vector<typename Dimension::Scalar> >& field) const {
  const vector<int>& gNodes = this->ghostNodes(field.nodeList());
  const unsigned nvoid = gNodes.size();
  for (unsigned k = 0; k < nvoid; ++k) field(gNodes[k]) = std::vector<Scalar>();
}

//------------------------------------------------------------------------------
// Provide the setViolationNodes for a NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
setViolationNodes(NodeList<Dimension>& nodeList) {
  this->addNodeList(nodeList);
}

//------------------------------------------------------------------------------
// Provide the updateViolationNodes for a NodeList, correcting positions and 
// H's.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
updateViolationNodes(NodeList<Dimension>& nodeList) {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for int fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
enforceBoundary(Field<Dimension, int>& field) const {
}

// Specialization for scalar fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Scalar>& field) const {
}

// Specialization for Vector fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
}

// Specialization for Tensor fields.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
}

// Specialization for third rank tensors.
template<typename Dimension>
void
CRKSPHVoidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
}

}
}

