//---------------------------------Spheral++----------------------------------//
// SPHHydroBaseRZ -- The SPH/ASPH hydrodynamic package for Spheral++,
//                   specialized for 2D RZ (cylindrical) geometry.
//
// Based on the methodology described in the dissertation
// de Catalunya Departament de Física i Enginyeria Nuclear, U. P., García Senz, D., (2012).
// AxisSPH:devising and validating an axisymmetric smoothed particle hydrodynamics code.
// TDX (Tesis Doctorals en Xarxa). Universitat Politècnica de Catalunya.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Tue Apr 26 16:06:21 PDT 2016
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "SPHHydroBaseRZ.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computeSPHOmegaGradhCorrection.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "FileIO/FileIO.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"

namespace Spheral {
namespace SPHSpace {

using namespace std;
using PhysicsSpace::GenericHydro;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FileIOSpace::FileIO;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using MeshSpace::Mesh;

using PhysicsSpace::MassDensityType;
using PhysicsSpace::HEvolutionType;

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SPHHydroBaseRZ::
SPHHydroBaseRZ(const SmoothingScaleBase<Dim<2> >& smoothingScaleMethod,
               ArtificialViscosity<Dim<2> >& Q,
               const TableKernel<Dim<2> >& W,
               const TableKernel<Dim<2> >& WPi,
               const double filter,
               const double cfl,
               const bool useVelocityMagnitudeForDt,
               const bool compatibleEnergyEvolution,
               const bool evolveTotalEnergy,
               const bool gradhCorrection,
               const bool XSPH,
               const bool correctVelocityGradient,
               const bool sumMassDensityOverAllNodeLists,
               const MassDensityType densityUpdate,
               const HEvolutionType HUpdate,
               const double epsTensile,
               const double nTensile,
               const Vector& xmin,
               const Vector& xmax):
  SPHHydroBase<Dim<2> >(smoothingScaleMethod,
                        Q,
                        W,
                        WPi,
                        filter,
                        cfl,
                        useVelocityMagnitudeForDt,
                        compatibleEnergyEvolution,
                        evolveTotalEnergy,
                        gradhCorrection,
                        XSPH,
                        correctVelocityGradient,
                        sumMassDensityOverAllNodeLists,
                        densityUpdate,
                        HUpdate,
                        epsTensile,
                        nTensile,
                        xmin,
                        xmax) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SPHHydroBaseRZ::
~SPHHydroBaseRZ() {
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SPHHydroBaseRZ::
evaluateDerivatives(const Dim<2>::Scalar time,
                    const Dim<2>::Scalar dt,
                    const DataBase<Dim<2> >& dataBase,
                    const State<Dim<2> >& state,
                    StateDerivatives<Dim<2> >& derivatives) const {

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  typedef Timing::Time Time;
  const double tiny = 1.0e-30;
  const Scalar W0 = W(0.0, 1.0);

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> omega = state.fields(HydroFieldNames::omegaGradh, 0.0);


  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  FieldList<Dimension, Tensor> localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  FieldList<Dimension, SymTensor> DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, Scalar> effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  FieldList<Dimension, Scalar> viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Scalar> XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const NodeList<Dimension>& nodeList = **itr;
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar hmin = nodeList.hmin();
    const Scalar hmax = nodeList.hmax();
    const Scalar hminratio = nodeList.hminratio();
    const int maxNumNeighbors = nodeList.maxNumNeighbors();
    const Scalar nPerh = nodeList.nodesPerSmoothingScale();

    // The scale for the tensile correction.
    const Scalar WnPerh = W(1.0/nPerh, 1.0);

    // Get the work field for this NodeList.
    Field<Dimension, Scalar>& workFieldi = nodeList.work();

    // Iterate over the internal nodes in this NodeList.
    for (ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Prepare to accumulate the time.
      const Time start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar& mi = mass(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const Scalar& epsi = specificThermalEnergy(nodeListi, i);
      const Scalar& Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar& ci = soundSpeed(nodeListi, i);
      const Scalar& omegai = omega(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar safeOmegai = safeInv(omegai, tiny);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      // RZ correction factors for node i.
      const Scalar zetai = (Hi*ri).y();
      const Scalar hrInvi = zetai*safeInvVar(ri.y());
      Scalar f1i, f2i, gradf1i, gradf2i;
      W.f1Andf2(zetai, f1i, f2i, gradf1i, gradf2i);
      gradf1i *= hrInvi;
      gradf2i *= hrInvi;
      const Scalar circi = 2.0*M_PI*ri.y();
      const Scalar circInvi = safeInvVar(circi);
      const Scalar rhoRZi = circi*rhoi;

      Scalar& rhoSumi = rhoSum(nodeListi, i);
      Scalar& normi = normalization(nodeListi, i);
      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      Tensor& Mi = M(nodeListi, i);
      Tensor& localMi = localM(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      Scalar& effViscousPressurei = effViscousPressure(nodeListi, i);
      Scalar& viscousWorki = viscousWork(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Scalar& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const double fweightij = 1.0; // (nodeListi == nodeListj ? 1.0 : 0.2);
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j
              const Vector& rj = position(nodeListj, j);
              const Scalar& mj = mass(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar& rhoj = massDensity(nodeListj, j);
              const Scalar& epsj = specificThermalEnergy(nodeListj, j);
              const Scalar& Pj = pressure(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar& cj = soundSpeed(nodeListj, j);
              const Scalar& omegaj = omega(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar safeOmegaj = safeInv(omegaj, tiny);
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Hdetj > 0.0);

              // RZ correction factors for node j.
              const Scalar zetaj = (Hj*rj).y();
              const Scalar hrInvj = zetaj*safeInvVar(rj.y());
              Scalar f1j, f2j, gradf1j, gradf2j;
              W.f1Andf2(zetaj, f1j, f2j, gradf1j, gradf2j);
              gradf1j *= hrInvj;
              gradf2j *= hrInvj;
              const Scalar circj = 2.0*M_PI*rj.y();
              const Scalar circInvj = safeInvVar(circj);
              const Scalar rhoRZj = circj*rhoj;

              Scalar& rhoSumj = rhoSum(nodeListj, j);
              Scalar& normj = normalization(nodeListj, j);
              Vector& DxDtj = DxDt(nodeListj, j);
              Scalar& DrhoDtj = DrhoDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Tensor& Mj = M(nodeListj, j);
              Tensor& localMj = localM(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              Scalar& effViscousPressurej = effViscousPressure(nodeListj, j);
              Scalar& viscousWorkj = viscousWork(nodeListj, j);
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
              Scalar& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
              Vector& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              Scalar& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              SymTensor& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Node displacement.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);

              // Symmetrized kernel weight and gradient.
              const Vector Hetai = Hi*etai.unitVector();
              const std::pair<double, double> WWi = W.kernelAndGradValue(etaMagi, Hdeti);
              const Scalar Wi = WWi.first;
              const Scalar gWi = WWi.second;
              const Vector gradWi = gWi*Hetai;
              const Vector gradWQi = WQ.gradValue(etaMagi, Hdeti) * Hetai;

              const Vector Hetaj = Hj*etaj.unitVector();
              const std::pair<double, double> WWj = W.kernelAndGradValue(etaMagj, Hdetj);
              const Scalar Wj = WWj.first;
              const Scalar gWj = WWj.second;
              const Vector gradWj = gWj*Hetaj;
              const Vector gradWQj = WQ.gradValue(etaMagj, Hdetj) * Hetaj;

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double rij2 = rij.magnitude2();
              const SymTensor thpt = rij.selfdyad()/max(tiny, rij2*FastMath::square(Dimension::pownu12(rij2)));
              weightedNeighborSumi += fweightij*std::abs(gWi);
              weightedNeighborSumj += fweightij*std::abs(gWj);
              massSecondMomenti += fweightij*gradWi.magnitude2()*thpt;
              massSecondMomentj += fweightij*gradWj.magnitude2()*thpt;

              // Contribution to the sum density.
              if (nodeListi == nodeListj) {
                rhoSumi += mj*Wi;
                rhoSumj += mi*Wj;
                normi += mi/rhoRZi*Wi;
                normj += mj/rhoRZj*Wj;
              }

              // Mass density evolution.  Here we sum up DrhoZR/Dt actually, to be converted at the end.
              DrhoDti += mj*((f1i*vi.y() - f2i*vj.y())*gradWi.y() + (gradf1i*vi.y() - gradf2i*vj.y())*Wi + f1i*(vi.z() - vj.z())*gradWi.z());
              DrhoDtj += mi*((f1j*vj.y() - f2j*vi.y())*gradWj.y() + (gradf1j*vj.y() - gradf2j*vi.y())*Wj + f1j*(vj.z() - vi.z())*gradWj.z());

              // Compute the pair-wise artificial viscosity.
              const Vector vij = vi - vj;
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Vector Qacci = 0.5*(QPiij.first *gradWQi);
              const Vector Qaccj = 0.5*(QPiij.second*gradWQj);
              // const Scalar workQi = 0.5*(QPiij.first *vij).dot(gradWQi);
              // const Scalar workQj = 0.5*(QPiij.second*vij).dot(gradWQj);
              const Scalar workQi = vij.dot(Qacci);
              const Scalar workQj = vij.dot(Qaccj);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);
              effViscousPressurei += mj/rhoj * Qi * Wi;
              effViscousPressurej += mi/rhoi * Qj * Wj;
              viscousWorki += mj*workQi;
              viscousWorkj += mi*workQj;

              // Determine an effective pressure including a term to fight the tensile instability.
//             const Scalar fij = epsTensile*pow(Wi/(Hdeti*WnPerh), nTensile);
              const Scalar fij = mEpsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
              const Scalar Ri = fij*(Pi < 0.0 ? -Pi : 0.0);
              const Scalar Rj = fij*(Pj < 0.0 ? -Pj : 0.0);
              const Scalar Peffi = Pi + Ri;
              const Scalar Peffj = Pj + Rj;

              // Acceleration.
              CHECK(rhoRZi > 0.0);
              CHECK(rhoRZj > 0.0);
              const double Prhoi = safeOmegai*Peffi*ri.y()/(rhoRZi*rhoRZi);
              const double Prhoj = safeOmegaj*Peffj*rj.y()/(rhoRZj*rhoRZj);
              const Vector deltaDvDti = -mj*(2.0*M_PI*(Prhoi*f1i*gradWi + Prhoj*gradWj) + Qacci + Qaccj);
              const Vector deltaDvDtj = -mi*(2.0*M_PI*(Prhoj*f1j*gradWj + Prhoi*gradWi) + Qaccj + Qacci);
              DvDti += deltaDvDti;
              DvDtj += deltaDvDtj;

              // Specific thermal energy evolution.
              DepsDti += mj*((2.0*M_PI*Peffi*ri.y()/(rhoRZi*rhoRZi)*
                              (f1i*vi.y() - f2i*vj.y())*gradWi.y() + f1i*(vi.x() - vj.x())*gradWi.x() + (gradf1i*vi.y() - gradf2i*vj.y())*Wi) + workQi);
              DepsDtj += mi*((2.0*M_PI*Peffj*rj.y()/(rhoRZj*rhoRZj)*
                              (f1j*vj.y() - f2j*vi.y())*gradWj.y() + f1j*(vj.x() - vi.x())*gradWj.x() + (gradf1j*vj.y() - gradf2j*vi.y())*Wj) + workQj);
              if (mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(deltaDvDti);
                pairAccelerationsj.push_back(deltaDvDtj);
              }

              // Velocity gradient.
              const Tensor deltaDvDxi = mj*vij.dyad(gradWi);
              const Tensor deltaDvDxj = mi*vij.dyad(gradWj);
              DvDxi -= deltaDvDxi;
              DvDxj -= deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi -= deltaDvDxi;
                localDvDxj -= deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (mXSPH and (nodeListi == nodeListj)) {
                const double fXSPH = max(0.0, min(1.0, abs(vij.dot(rij)*safeInv(vij.magnitude()*rij.magnitude()))));
                CHECK(fXSPH >= 0.0 and fXSPH <= 1.0);
                XSPHWeightSumi += fXSPH*mj/rhoRZj*Wi;
                XSPHWeightSumj += fXSPH*mi/rhoRZi*Wj;
                XSPHDeltaVi -= fXSPH*mj/rhoRZj*Wi*vij;
                XSPHDeltaVj += fXSPH*mi/rhoRZi*Wj*vij;
              }

              // Linear gradient correction term.
              Mi -= mj*rij.dyad(gradWi);
              Mj -= mi*rij.dyad(gradWj);
              if (nodeListi == nodeListj) {
                localMi -= mj*rij.dyad(gradWi);
                localMj -= mi*rij.dyad(gradWj);
              }
            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not mCompatibleEnergyEvolution or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Get the time for pairwise interactions.
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // Add the self-contribution to density sum.
      rhoSumi = (rhoSumi + mi*W0*Hdeti)*f1i*circInvi;
      normi += mi/rhoRZi*W0*Hdeti;

      // Finish the acceleration, adding the hoop terms.
      DvDti.y(DvDti.y() + 2.0*M_PI*Pi/rhoRZi*(1.0 - ri.y()*safeInvVar(f1i)*gradf1i));

      // Finish the gradient of the velocity.
      CHECK(rhoRZi > 0.0);
      if (this->mCorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoRZi;
      }
      if (this->mCorrectVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      } else {
        localDvDxi /= rhoRZi;
      }

      // Finish the continuity equation.
      DrhoDti *= (ri.y() - DvDti.y())*safeInvVar(2.0*M_PI*FastMath::square(ri.y()));

      // If needed finish the total energy derivative.
      if (mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        XSPHWeightSumi += Hdeti*mi/rhoRZi*W0;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDti = vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi);
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             ri,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                        ri,
                                                        weightedNeighborSumi,
                                                        massSecondMomenti,
                                                        W,
                                                        hmin,
                                                        hmax,
                                                        hminratio,
                                                        nPerh,
                                                        connectivityMap,
                                                        nodeListi,
                                                        i);

      // Increment the work for i.
      worki += Timing::difference(start, Timing::currentTime());

      // Now add the pairwise time for each neighbor we computed here.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
          Field<Dimension, Scalar>& workFieldj = nodeLists[nodeListj]->work();
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              workFieldj(j) += deltaTimePair;
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
void
SPHHydroBaseRZ::
finalize(const Dim<2>::Scalar time,
         const Dim<2>::Scalar dt,
         DataBase<Dim<2> >& dataBase,
         State<Dim<2> >& state,
         StateDerivatives<Dim<2> >& derivs) {

  // Base class finalization does most of the work.
  SPHHydroBase<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // If we're doing the rigorous sum, we need to apply RZ correction factor and
  // convert to 3D density.
  if (densityUpdate() == PhysicsSpace::RigorousSumDensity or
      densityUpdate() == PhysicsSpace::CorrectedSumDensity) {
    const TableKernel<Dimension>& W = this->kernel();
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const unsigned numNodeLists = massDensity.numFields();
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = massDensity[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        const Vector& xi = position(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);
        const Scalar zetai = (Hi*xi).y();
        const Scalar fi = W.f1(zetai);
        massDensity(nodeListi, i) *= fi*safeInvVar(2.0*M_PI*xi.y());
      }
    }
  }
}

}
}

