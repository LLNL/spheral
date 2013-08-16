//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructred the 
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#include "FiniteVolumeViscosity.hh"
#include "DataOutput/Restart.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Material/EquationOfState.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "Mesh/Mesh.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using DataOutput::Restart;
using FieldSpace::Field;
using FieldSpace::FieldList;
using DataBaseSpace::DataBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NeighborSpace::Neighbor;
using Material::EquationOfState;
using BoundarySpace::Boundary;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using MeshSpace::Mesh;

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
FiniteVolumeViscosity<Dimension>::
FiniteVolumeViscosity(const Scalar Clinear,
                      const Scalar Cquadratic,
                      const bool scalar):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic),
  mScalar(scalar),
  mDvDx(FieldList<Dimension, Tensor>::Copy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FiniteVolumeViscosity<Dimension>::
~FiniteVolumeViscosity() {
}

//------------------------------------------------------------------------------
// The required method to compute the artificial viscous P/rho^2.
//------------------------------------------------------------------------------
template<typename Dimension>
pair<typename Dimension::Tensor,
     typename Dimension::Tensor>
FiniteVolumeViscosity<Dimension>::
Piij(const unsigned nodeListi, const unsigned i, 
     const unsigned nodeListj, const unsigned j,
     const Vector& xi,
     const Vector& etai,
     const Vector& vi,
     const Scalar rhoi,
     const Scalar csi,
     const SymTensor& Hi,
     const Vector& xj,
     const Vector& etaj,
     const Vector& vj,
     const Scalar rhoj,
     const Scalar csj,
     const SymTensor& Hj) const {

  const double Cl = this->mClinear;
  const double Cq = this->mCquadratic;
  const double eps2 = this->mEpsilon2;
  const bool balsaraShearCorrection = this->mBalsaraShearCorrection;

  // Are we applying the shear corrections?
  const Vector vij = vi - vj;
  const Scalar fshear = (balsaraShearCorrection ? 
                         abs(vij.unitVector().dot(vij.unitVector())) :
                         1.0);

  const Vector xji = xj - xi;
  const Vector xjihat = xji.unitVector();
  const Scalar hi = 1.0/(Hi*xjihat).magnitude();
  const Scalar hj = 1.0/(Hj*xjihat).magnitude();
  const Scalar DvDxi = min(0.0, mDvDx(nodeListi, i).Trace());
  const Scalar DvDxj = min(0.0, mDvDx(nodeListj, j).Trace());
  const Scalar Pii = fshear*(-Cl*csi*DvDxi + Cq*hi*DvDxi*DvDxi)*hi/rhoi;
  const Scalar Pij = fshear*(-Cl*csj*DvDxj + Cq*hj*DvDxj*DvDxj)*hj/rhoj;
  return make_pair(Pii*Tensor::one, Pij*Tensor::one);
}

//------------------------------------------------------------------------------
// Initialize for the FluidNodeLists in the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FiniteVolumeViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           const State<Dimension>& state,
           const StateDerivatives<Dimension>& derivs,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
	   const typename Dimension::Scalar time,
	   const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

  typedef typename Mesh<Dimension>::Zone Zone;
  typedef typename Mesh<Dimension>::Face Face;

  // Call the ancestor.
  ArtificialViscosity<Dimension>::initialize(dataBase,
                                             state,
                                             derivs,
                                             boundaryBegin,
                                             boundaryEnd,
                                             time,
                                             dt,
                                             W);

  // Prepare our result.
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, "FV DvDx", true);

  // Make a finite-volume estimate of the local (to each Voronoi cell) velocity
  // gradient.
  unsigned nodeListj, j;
  Scalar Vi, Vj;
  const Mesh<Dimension>& mesh = state.mesh();
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const unsigned numNodeLists = velocity.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = velocity[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Vector& vi = velocity(nodeListi, i);
      Tensor& DvDxi = mDvDx(nodeListi, i);
      const Zone& zonei = mesh.zone(nodeListi, i);
      const vector<int>& faces = zonei.faceIDs();
      Vi = zonei.volume();
      for (vector<int>::const_iterator fitr = faces.begin();
           fitr != faces.end();
           ++fitr) {
        const Face& faceij = mesh.face(*fitr);
        const int oppZoneID = faceij.oppositeZoneID(zonei.ID());
        if (Mesh<Dimension>::positiveID(oppZoneID) == Mesh<Dimension>::UNSETID) {
          nodeListj = nodeListi;
          j = i;
        } else {
          mesh.lookupNodeListID(Mesh<Dimension>::positiveID(oppZoneID), nodeListj, j);
        }
        const Vector& vj = velocity(nodeListj, j);
        const Vector vij = 0.5*(vi + vj);
        const Vector dA = faceij.area() * faceij.unitNormal() * sgn(*fitr);
        DvDxi -= vij*dA;
      }
      DvDxi /= Vi;
    }
  }

  // // Apply boundary conditions to our intermediate estimate.
  // for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
  //      boundItr != boundaryEnd;
  //      ++boundItr) {
  //   (*boundItr)->applyFieldListGhostBoundary(mDvDx);
  // }

  // for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
  //      boundItr != boundaryEnd;
  //      ++boundItr) {
  //   (*boundItr)->finalizeGhostBoundary();
  // }

  // // Make a simple smoothing iteration over DvDx.
  // FieldList<Dimension, Tensor> DvDx_smooth(mDvDx);
  // DvDx_smooth.copyFields();
  // for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //   const unsigned n = velocity[nodeListi]->numInternalElements();
  //   for (unsigned i = 0; i != n; ++i) {
  //     const Zone& zonei = mesh.zone(nodeListi, i);
  //     const vector<int>& faces = zonei.faceIDs();
  //     Vi = zonei.volume();
  //     Scalar Vtot = Vi;
  //     DvDx_smooth(nodeListi, i) = Vi*mDvDx(nodeListi, i);
  //     for (vector<int>::const_iterator fitr = faces.begin();
  //          fitr != faces.end();
  //          ++fitr) {
  //       const Face& faceij = mesh.face(*fitr);
  //       const int oppZoneID = faceij.oppositeZoneID(zonei.ID());
  //       if (Mesh<Dimension>::positiveID(oppZoneID) == Mesh<Dimension>::UNSETID) {
  //         nodeListj = nodeListi;
  //         j = i;
  //         Vj = Vi;
  //       } else {
  //         mesh.lookupNodeListID(Mesh<Dimension>::positiveID(oppZoneID), nodeListj, j);
  //         Vj = mesh.zone(Mesh<Dimension>::positiveID(oppZoneID)).volume();
  //       }
  //       Vtot += Vj;
  //       DvDx_smooth(nodeListi, i) += Vj*mDvDx(nodeListj, j);
  //     }
  //     DvDx_smooth(nodeListi, i) /= Vtot;
  //   }
  // }
  // mDvDx = DvDx_smooth;

  // Apply boundary conditions.  We depend on someone else finalizing these 
  // between now and calls to Piij.
  for (typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundItr = boundaryBegin;
       boundItr != boundaryEnd;
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(mDvDx);
  }

}


//------------------------------------------------------------------------------
// Are we enforcing scalar Q?
//------------------------------------------------------------------------------
template<typename Dimension>
bool
FiniteVolumeViscosity<Dimension>::
scalar() const {
  return mScalar;
}

//------------------------------------------------------------------------------
// The finite-volume gradient of the velocity.
//------------------------------------------------------------------------------
template<typename Dimension>
const FieldList<Dimension, typename Dimension::Tensor>&
FiniteVolumeViscosity<Dimension>::
DvDx() const {
  return mDvDx;
}

}
}
