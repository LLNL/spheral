//---------------------------------Spheral++----------------------------------//
// An experimental hour glass control algorithm for SPH, based on the using
// center of mass estimates in Voronoi cells.
//
// Created by JMO, Tue Jun 28 14:54:03 PDT 2011
//----------------------------------------------------------------------------//
#include "Hydro/VoronoiHourglassControl.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Mesh/Mesh.hh"
#include "Mesh/MeshPolicy.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "CRKSPH/computeCRKSPHCorrections.hh"
#include "FieldOperations/monotonicallyLimitedGradient.hh"
#include "Distributed/Communicator.hh"
#include "Distributed/allReduce.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {


//------------------------------------------------------------------------------
// Find the center of mass for a cell given a slope for the density in the cell.
//------------------------------------------------------------------------------
// 1D
inline
Dim<1>::Vector
centerOfMass(const Dim<1>::Vector& xi,
             const Dim<1>::Scalar& rhoi,
             const Dim<1>::Vector& gradRhoi,
             const Mesh<Dim<1> >& mesh,
             const unsigned zoneID) {
  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  REQUIRE(zoneID < mesh.numZones());

  // Read the node positions.
  const vector<unsigned>& nodeIDs = mesh.zone(zoneID).nodeIDs();
  CHECK(nodeIDs.size() == 2);
  const Scalar x0 = mesh.node(nodeIDs[0]).position().x();
  const Scalar x1 = mesh.node(nodeIDs[1]).position().x();
  const Scalar xxi = xi.x();
  CHECK(x0 < x1);
  CHECK2(x0 <= xxi and xxi <= x1, zoneID << " : " << x0 << " " << xxi << " " << x1);

  const Scalar dx  = x1 - x0;
  const Scalar dx2 = x1*x1 - x0*x0;
  const Scalar dx3 = x1*x1*x1 - x0*x0*x0;
  const Scalar grhoi = gradRhoi.x();
  const Scalar A = rhoi - grhoi*xxi;
  const Scalar num = 0.5*A*dx2 + grhoi*dx3/3.0;
  const Scalar den = A*dx + 0.5*grhoi*dx2;
  CHECK(den > 0.0);
  // cerr << "  --> " << zoneID << " " << (num/den - xxi)/0.01 << endl;
  return Vector(num/den);
}

// 2D
inline
Dim<2>::Vector
centerOfMass(const Dim<2>::Vector& xi,
             const Dim<2>::Scalar& rhoi,
             const Dim<2>::Vector& gradRhoi,
             const Mesh<Dim<2> >& mesh,
             const unsigned zoneID) {
  typedef Dim<2>::Vector Vector;
  typedef Mesh<Dim<2> >::Zone Zone;
  const double onethird = 1.0/3.0;

  // Loop over the vertices, decompose into triangles, and accumulate
  // the weighted center of mass.

  // Extract the zone info.
  const Zone& zone = mesh.zone(zoneID);
  const vector<unsigned>& nodeIDs = zone.nodeIDs();
  const unsigned numNodes = nodeIDs.size();
  const Vector zpos = zone.position();

  // Walk the triangles and compute their weighted contribution to the 
  // overall center of mass.
  unsigned i, j;
  Vector result, tri, nj, ni = mesh.node(nodeIDs[0]).position();
  double weightSum = numeric_limits<double>::min(), weight;
  for (i = 0; i != numNodes; ++i) {
    j = (i + 1) % numNodes;
    nj = mesh.node(nodeIDs[j]).position();
    tri = onethird*(zpos + ni + nj);
    weight = max(0.0, (rhoi + (tri - xi).dot(gradRhoi)) * ((ni - zpos).cross(nj - zpos).z()));
    CHECK2(weight >= 0.0,
           "Bad weight: " << weight << " " << rhoi << " " << gradRhoi << " " 
           << (rhoi + (tri - xi).dot(gradRhoi)) << " "<< ((ni - zpos).cross(nj - zpos).z()));
    weightSum += weight;
    result += weight * tri;
    ni = nj;
  }
  CHECK2(weightSum > 0.0, weightSum << " " << rhoi << " " << gradRhoi << xi);
  result /= weightSum;

  // That's it.
  ENSURE2(zone.convexHull().convexContains(result),
          zpos << " " << xi << " " << result << " : " << zone.convexHull().distance(result));
  return result;
}

// 3D
inline
Dim<3>::Vector
centerOfMass(const Dim<3>::Vector& xi,
             const Dim<3>::Scalar& rhoi,
             const Dim<3>::Vector& gradRhoi,
             const Mesh<Dim<3> >& mesh,
             const unsigned zoneID) {
  VERIFY(false);
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiHourglassControl<Dimension>::
VoronoiHourglassControl(const TableKernel<Dimension>& W,
                        const unsigned order,
                        const unsigned limiter,
                        const double fraction,
                        const FieldList<Dimension, int>& mask):
  Physics<Dimension>(),
  mW(W),
  mOrder(order),
  mLimiter(limiter),
  mFraction(fraction),
  mMask(mask),
  mA(FieldStorageType::CopyFields),
  mWeight(FieldStorageType::CopyFields),
  mGradRho(FieldStorageType::CopyFields),
  mB(FieldStorageType::CopyFields),
  mC(FieldStorageType::CopyFields),
  mGradA(FieldStorageType::CopyFields),
  mD(FieldStorageType::CopyFields),
  mGradB(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
VoronoiHourglassControl<Dimension>::
~VoronoiHourglassControl() {
}

//------------------------------------------------------------------------------
// Determine the principle derivatives for the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiHourglassControl<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename VoronoiHourglassControl<Dimension>::TimeStepType
VoronoiHourglassControl<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const typename Dimension::Scalar currentTime) const {
  return TimeStepType(FLT_MAX, "No vote.");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiHourglassControl<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Register the Voronoi mesh.
  PolicyPointer meshPolicy(new MeshPolicy<Dimension>(*this));
  state.enroll(HydroFieldNames::mesh, meshPolicy);

  // These are nice to register for analysis.
  dataBase.resizeFluidFieldList(mMask, 1, HydroFieldNames::hourglassMask, false);
  dataBase.resizeFluidFieldList(mGradRho, Vector::zero, "grad " + HydroFieldNames::massDensity);
  state.enroll(mMask);
  state.enroll(mGradRho);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiHourglassControl<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Finalize -- filter the positions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
VoronoiHourglassControl<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Extract our state.
  const Mesh<Dimension>& mesh = state.mesh();
  FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Scalar>& mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Scalar>& rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, SymTensor>& H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Tensor>& gradv = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);

  // Flatten the density to a 1D field corresponding to the mesh zones.
  const unsigned numNodeLists = position.size();
  unsigned nodeListi, nodeListj, i, j, k, kk, ii, offset, firstGhostNode;
  Scalar Ai, Wj, gWj, rhoMin = numeric_limits<Scalar>::max(), rhoMax = -1.0;
  Vector ri, rij, etaj, Bi, gradAi, gradWj, rNode, thpt;
  Tensor gradBi;
  SymTensor Hj;
  vector<Vector> rNodes;
  vector<Scalar> rhoNodes;
  vector<Scalar> rhoZones;
  rhoZones.reserve(mesh.numZones());
  for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (i = 0; i != rho[nodeListi]->numInternalElements(); ++i) {
      rhoZones.push_back(rho(nodeListi, i));
      rhoMin = std::min(rhoMin, rhoZones.back());
      rhoMax = std::max(rhoMax, rhoZones.back());
    }
  }
  CHECK(rhoZones.size() == mesh.numZones());
  rhoMin = allReduce(rhoMin, SPHERAL_OP_MIN);
  rhoMax = allReduce(rhoMax, SPHERAL_OP_MAX);

  // Compute the CRKSPH limited gradient of the density if we're doing first order.
  if (mOrder > 0) {
    dataBase.resizeFluidFieldList(mA, 0.0,              HydroFieldNames::A_CRKSPH);
    dataBase.resizeFluidFieldList(mB, Vector::zero,     HydroFieldNames::B_CRKSPH);
    dataBase.resizeFluidFieldList(mC, Vector::zero,     HydroFieldNames::C_CRKSPH);
    dataBase.resizeFluidFieldList(mD, Tensor::zero,     HydroFieldNames::D_CRKSPH);
    dataBase.resizeFluidFieldList(mGradA, Vector::zero, HydroFieldNames::gradA_CRKSPH);
    dataBase.resizeFluidFieldList(mGradB, Tensor::zero, HydroFieldNames::gradB_CRKSPH);
    mWeight = mass/rho;
    for (i = 0; i != mA.size(); ++i) {
      state.enroll(*mA[i]);
      state.enroll(*mB[i]);
      state.enroll(*mC[i]);
      state.enroll(*mD[i]);
      state.enroll(*mGradA[i]);
      state.enroll(*mGradB[i]);
      state.enroll(*mWeight[i]);
    }

    // Compute the CRKSPH correction terms.
    dataBase.updateConnectivityMap();
    const ConnectivityMap<Dimension>& cm = dataBase.connectivityMap();
    computeCRKSPHCorrections(cm, mW, mWeight, position, H, 
                                      mA, mA, mB, mC, mD, mGradA, mGradB);

    // Find the gradient of the density.
    for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      firstGhostNode = cm.nodeLists()[nodeListi]->firstGhostNode();
      for (typename ConnectivityMap<Dimension>::const_iterator iItr = cm.begin(nodeListi);
           iItr != cm.end(nodeListi);
           ++iItr) {
        i = *iItr;
        offset = mesh.offset(nodeListi);
        ri = position(nodeListi, i);
        Ai = mA(nodeListi, i);
        Bi = mB(nodeListi, i);
        gradAi = mGradA(nodeListi, i);
        gradBi = mGradB(nodeListi, i);
        const vector<vector<int> >& fullConnectivity = cm.connectivityForNode(nodeListi, i);
        CHECK(fullConnectivity.size() == numNodeLists);
        for (nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          const vector<int>& connectivity = fullConnectivity[nodeListj];
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            j = *jItr;
            Hj = H(nodeListj, j);
            rij = ri - position(nodeListj, j);
            etaj = Hj*rij;
            CRKSPHKernelAndGradient(mW, rij, etaj, Hj, Hj.Determinant(), 
                                             Ai, Bi, gradAi, gradBi,
                                             Wj, gWj, gradWj);
            mGradRho(nodeListi, i) += rho(nodeListj, j) * mWeight(nodeListj, j)*gradWj;
          }
        }

        if (mLimiter > 0 and i < firstGhostNode) {
          // Apply monotonic limiting to the gradient.
          rNodes = vector<Vector>();
          rhoNodes = vector<Scalar>();
          CHECK2(offset + i < mesh.numZones(), "Bad zone indexing:  " << offset << " " << i << " " << mesh.numZones());
          const vector<unsigned>& nodeIDs = mesh.zone(offset + i).nodeIDs();
          for (k = 0; k != nodeIDs.size(); ++k) {
            j = nodeIDs[k];
            const vector<unsigned>& zoneIDs = mesh.node(j).zoneIDs();
            rNode = mesh.node(j).position();
            for (kk = 0; kk != zoneIDs.size(); ++kk) {
              ii = zoneIDs[kk];
              if (ii != offset + i and ii != Mesh<Dimension>::UNSETID) {
                CHECK2(ii < rhoZones.size(), 
                       "Bad zone ID:  " << ii << " " << mesh.numZones() << " " << rhoZones.size() << " " << Mesh<Dimension>::UNSETID << endl
                       << ri << " " << rNode);
                rNodes.push_back(rNode);
                rhoNodes.push_back(rhoZones[ii]);
              }
            }
          }
          if (mLimiter == 1) {
            mGradRho(nodeListi, i) = scalarLimitedGradient<Dimension, Scalar>(rho(nodeListi, i), 
                                                                                          mGradRho(nodeListi, i),
                                                                                          ri, rNodes, rhoNodes);
          } else {
            CHECK(mLimiter == 2);
            mGradRho(nodeListi, i) = tensorLimitedGradient<Dimension, Scalar>(rho(nodeListi, i), 
                                                                                          mGradRho(nodeListi, i),
                                                                                          ri, rNodes, rhoNodes);
          }
        }
      }
    }
  }

  // Walk the NodeLists.
  offset = 0;
  const Scalar drhoInv = safeInv(rhoMax - rhoMin);
  Scalar dv;
  Vector deltacom;
  for (unsigned nodeListi = 0; nodeListi != position.size(); ++nodeListi) {
    const unsigned numNodes = position[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != numNodes; ++i) {
      if (mMask(nodeListi, i) == 1) {
        Vector& xi = position(nodeListi, i);
        deltacom = centerOfMass(position(nodeListi, i),
                                rho(nodeListi, i),
                                mGradRho(nodeListi, i),
                                mesh,
                                offset + i) - xi;
        if (mFraction > 0.0) {
          dv = (gradv(nodeListi, i)*deltacom).magnitude();
          xi += min(deltacom.magnitude(), mFraction*dv*dt)*(deltacom.unitVector());
        } else {
          xi += deltacom;
        }
        // const double f = 0.5;
        // position(nodeListi, i) = (f*position(nodeListi, i) + (1.0 - f)*centerOfMass(position(nodeListi, i),
        //                                                                             rho(nodeListi, i),
        //                                                                             mGradRho(nodeListi, i),
        //                                                                             mesh,
        //                                                                             offset + i));
      }
    }
    offset += numNodes;
  }
}

}
