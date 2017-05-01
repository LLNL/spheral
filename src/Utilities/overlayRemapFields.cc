//------------------------------------------------------------------------------
// Use geometric clipping to remap a set of conserved fields.
// Currently only works single NodeList -> single NodeList, no boundaries.
//------------------------------------------------------------------------------

#include <map>
#include <algorithm>

#include "overlayRemapFields.hh"
#include "r3d_utils.hh"
#include "DataBase/DataBase.hh"
#include "CRKSPH/computeVoronoiVolume.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using DataBaseSpace::DataBase;
using NeighborSpace::ConnectivityMap;
using BoundarySpace::Boundary;
using NeighborSpace::Neighbor;

template<typename Dimension>
void
overlayRemapFields(const vector<Boundary<Dimension>*>& boundaries,
                   const vector<Field<Dimension, typename Dimension::Scalar>*>& scalarDonorFields,
                   const vector<Field<Dimension, typename Dimension::Vector>*>& vectorDonorFields,
                   const vector<Field<Dimension, typename Dimension::Tensor>*>& tensorDonorFields,
                   const vector<Field<Dimension, typename Dimension::SymTensor>*>& symTensorDonorFields,
                   vector<Field<Dimension, typename Dimension::Scalar>*>& scalarAcceptorFields,
                   vector<Field<Dimension, typename Dimension::Vector>*>& vectorAcceptorFields,
                   vector<Field<Dimension, typename Dimension::Tensor>*>& tensorAcceptorFields,
                   vector<Field<Dimension, typename Dimension::SymTensor>*>& symTensorAcceptorFields) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  typedef typename FacetedVolume::Facet Facet;
  typedef GeomPlane<Dimension> Plane;

  // Preconditions.
  VERIFY2(scalarDonorFields.size() == scalarAcceptorFields.size(), "overlayRemapFields ERROR: number of acceptor scalar fields does not match number of donors.");
  VERIFY2(vectorDonorFields.size() == vectorAcceptorFields.size(), "overlayRemapFields ERROR: number of acceptor vector fields does not match number of donors.");
  VERIFY2(tensorDonorFields.size() == tensorAcceptorFields.size(), "overlayRemapFields ERROR: number of acceptor tensor fields does not match number of donors.");
  VERIFY2(symTensorDonorFields.size() == symTensorAcceptorFields.size(), "overlayRemapFields ERROR: number of acceptor symTensor fields does not match number of donors.");
  const unsigned nScalarFields = scalarDonorFields.size();
  const unsigned nVectorFields = vectorDonorFields.size();
  const unsigned nTensorFields = tensorDonorFields.size();
  const unsigned nSymTensorFields = symTensorDonorFields.size();

  // Find the donor and acceptor NodeLists.
  const NodeList<Dimension> *donorNodeListPtr = NULL, *acceptorNodeListPtr = NULL;
  for (const Field<Dimension, Scalar>* field: scalarDonorFields) {
    VERIFY2(donorNodeListPtr == NULL or donorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all donor fields must be on same NodeList.");
    donorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Scalar>* field: scalarDonorFields) {
    VERIFY2(donorNodeListPtr == NULL or donorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all donor fields must be on same NodeList.");
    donorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Scalar>* field: scalarDonorFields) {
    VERIFY2(donorNodeListPtr == NULL or donorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all donor fields must be on same NodeList.");
    donorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Scalar>* field: scalarDonorFields) {
    VERIFY2(donorNodeListPtr == NULL or donorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all donor fields must be on same NodeList.");
    donorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Scalar>* field: scalarAcceptorFields) {
    VERIFY2(acceptorNodeListPtr == NULL or acceptorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all acceptor fields must be on same NodeList.");
    acceptorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Scalar>* field: scalarAcceptorFields) {
    VERIFY2(acceptorNodeListPtr == NULL or acceptorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all acceptor fields must be on same NodeList.");
    acceptorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Scalar>* field: scalarAcceptorFields) {
    VERIFY2(acceptorNodeListPtr == NULL or acceptorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all acceptor fields must be on same NodeList.");
    acceptorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Scalar>* field: scalarAcceptorFields) {
    VERIFY2(acceptorNodeListPtr == NULL or acceptorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all acceptor fields must be on same NodeList.");
    acceptorNodeListPtr = field->nodeListPtr();
  }
  Neighbor<Dimension>& neighborD = donorNodeListPtr->neighbor();
  Neighbor<Dimension>& neighborA = acceptorNodeListPtr->neighbor();

  // Build the donor volumes.
  Field<Dimension, FacetedVolume> donorCells("donor cells", *donorNodeListPtr);
  {
    DataBase<Dimension> db;
    db.appendNodeList(*dynamic_cast<FluidNodeList<Dimension>*>(const_cast<NodeList<Dimension>*>(donorNodeListPtr)));
    neighborD.updateNodes();
    for (Boundary<Dimension>* boundPtr: boundaries) boundPtr->setAllGhostNodes(db);
    for (Boundary<Dimension>* boundPtr: boundaries) boundPtr->finalizeGhostBoundary();
    db.updateConnectivityMap(false);
    const ConnectivityMap<Dimension>& cm = db.connectivityMap();
    const FieldList<Dimension, Vector> position = db.fluidPosition();
    const FieldList<Dimension, SymTensor> H = db.fluidHfield();
    const FieldList<Dimension, Scalar> rho = db.fluidMassDensity();
    const FieldList<Dimension, Vector> gradrho = db.newFluidFieldList(Vector::zero, "rho gradient");
    const FieldList<Dimension, Scalar> weight = db.newFluidFieldList(1.0, "weight");
    FieldList<Dimension, int> surfacePoint = db.newFluidFieldList(0, "surface point");
    FieldList<Dimension, Scalar> vol = db.newFluidFieldList(0.0, "volume");
    FieldList<Dimension, Vector> deltaMedian = db.newFluidFieldList(Vector::zero, "displacement");
    FieldList<Dimension, FacetedVolume> cells_fl(FieldSpace::FieldStorageType::Reference);
    cells_fl.appendField(donorCells);
    CRKSPHSpace::computeVoronoiVolume(position, H, rho, gradrho, cm, 2.0, vector<FacetedVolume>(), vector<vector<FacetedVolume>>(), weight,
                                      surfacePoint, vol, deltaMedian, cells_fl);
  }
    
  // Build the acceptor volumes.
  Field<Dimension, FacetedVolume> acceptorCells("acceptor cells", *acceptorNodeListPtr);
  {
    neighborA.updateNodes();
    DataBase<Dimension> db;
    db.appendNodeList(*dynamic_cast<FluidNodeList<Dimension>*>(const_cast<NodeList<Dimension>*>(acceptorNodeListPtr)));
    db.updateConnectivityMap(false);
    const ConnectivityMap<Dimension>& cm = db.connectivityMap();
    const FieldList<Dimension, Vector> position = db.fluidPosition();
    const FieldList<Dimension, SymTensor> H = db.fluidHfield();
    const FieldList<Dimension, Scalar> rho = db.fluidMassDensity();
    const FieldList<Dimension, Vector> gradrho = db.newFluidFieldList(Vector::zero, "rho gradient");
    const FieldList<Dimension, Scalar> weight = db.newFluidFieldList(1.0, "weight");
    FieldList<Dimension, int> surfacePoint = db.newFluidFieldList(0, "surface point");
    FieldList<Dimension, Scalar> vol = db.newFluidFieldList(0.0, "volume");
    FieldList<Dimension, Vector> deltaMedian = db.newFluidFieldList(Vector::zero, "displacement");
    FieldList<Dimension, FacetedVolume> cells_fl(FieldSpace::FieldStorageType::Reference);
    cells_fl.appendField(acceptorCells);
    CRKSPHSpace::computeVoronoiVolume(position, H, rho, gradrho, cm, 2.0, vector<FacetedVolume>(), vector<vector<FacetedVolume>>(), weight,
                                      surfacePoint, vol, deltaMedian, cells_fl);
  }

  // Look for intersecting node volumes.
  const unsigned nD = donorNodeListPtr->numInternalNodes(), nA = acceptorNodeListPtr->numInternalNodes();
  const Field<Dimension, Vector>& posD = donorNodeListPtr->positions();
  const Field<Dimension, SymTensor>& HD = donorNodeListPtr->Hfield();
  const Field<Dimension, Vector>& posA = acceptorNodeListPtr->positions();
  const Field<Dimension, SymTensor>& HA = acceptorNodeListPtr->Hfield();
  Field<Dimension, vector<unsigned>> intersectIndices("intersection indices", *donorNodeListPtr);
  Field<Dimension, vector<Scalar>> intersectVols("intesection volumes", *donorNodeListPtr);
  for (unsigned i = 0; i != nD; ++i) {
    neighborA.setMasterList(posD(i), HD(i));
    neighborA.setRefineNeighborList(posD(i), HD(i));
    for (typename Neighbor<Dimension>::const_iterator jitr = neighborA.refineNeighborBegin();
         jitr != neighborA.refineNeighborEnd();
         ++jitr) {
      const int j = *jitr;
      if (donorCells(i).intersect(acceptorCells(j))) {
        const vector<Facet>& facets = acceptorCells(j).facets();
        vector<Plane> planes;
        planes.reserve(facets.size());
        for (const Facet& facet: facets) planes.push_back(Plane(facet.position(), -facet.normal()));
        const Scalar Vi = clipFacetedVolume(donorCells(i), planes).volume();
        if (Vi > 0.0) {
          intersectIndices(i).push_back(j);
          intersectVols(i).push_back(Vi);
        }
      }
    }
  }

  // Now we can go through and splat the conserved values from the donor to acceptor volumes.
  for (unsigned i = 0; i != nD; ++i) {
    const unsigned n = intersectIndices(i).size();
    CHECK(intersectVols(i).size() == n);
    const Scalar voltotInv = safeInv(accumulate(intersectVols(i).begin(), intersectVols(i).end(), 0.0));
    for (unsigned k = 0; k != n; ++k) {
      const unsigned j = intersectIndices(i)[k];
      const Scalar f = intersectVols(i)[k]*voltotInv;
      for (unsigned kk = 0; kk != nScalarFields; ++kk) (*scalarAcceptorFields[kk])(j) += f*(*scalarDonorFields[kk])(i);
      for (unsigned kk = 0; kk != nVectorFields; ++kk) (*vectorAcceptorFields[kk])(j) += f*(*vectorDonorFields[kk])(i);
      for (unsigned kk = 0; kk != nTensorFields; ++kk) (*tensorAcceptorFields[kk])(j) += f*(*tensorDonorFields[kk])(i);
      for (unsigned kk = 0; kk != nSymTensorFields; ++kk) (*symTensorAcceptorFields[kk])(j) += f*(*symTensorDonorFields[kk])(i);
    }
  }
}

}
