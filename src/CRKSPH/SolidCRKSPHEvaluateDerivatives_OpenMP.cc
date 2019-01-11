namespace Spheral {

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();
  const auto& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto XSPH = this->XSPH();
  const auto epsTensile = this->epsilonTensile();
  const auto order = this->correctionOrder();
  const auto correctionMin = this->correctionMin();
  const auto correctionMax = this->correctionMax();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto entropy = state.fields(HydroFieldNames::entropy, Scalar());
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  const auto gradDamage = state.fields(SolidFieldNames::damageGradient, Vector::zero);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  const auto B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  const auto C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  const auto gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  const auto gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  const auto gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);
  const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  const auto voidPoint = state.fields(HydroFieldNames::voidPoint, 0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(entropy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(gradDamage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(A.size() == numNodeLists);
  CHECK(B.size() == numNodeLists);
  CHECK(C.size() == numNodeLists or order != CRKOrder::QuadraticOrder);
  CHECK(gradA.size() == numNodeLists);
  CHECK(gradB.size() == numNodeLists);
  CHECK(gradC.size() == numNodeLists or order != CRKOrder::QuadraticOrder);
  CHECK(surfacePoint.size() == numNodeLists);
  CHECK(voidPoint.size() == numNodeLists);

  // Derivative FieldLists.
  auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto DSDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto gradRho = derivatives.fields(HydroFieldNames::massDensityGradient, Vector::zero);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);
  CHECK(gradRho.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (compatibleEnergy) {
    auto nodeListi = 0;
    for (auto itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Some scratch variables.
  Scalar Ai, Aj;
  Vector gradAi, gradAj, forceij, forceji;
  Vector Bi = Vector::zero, Bj = Vector::zero;
  Tensor Ci = Tensor::zero, Cj = Tensor::zero;
  Tensor gradBi = Tensor::zero, gradBj = Tensor::zero;
  ThirdRankTensor gradCi = ThirdRankTensor::zero, gradCj = ThirdRankTensor::zero;
  Scalar gWi, gWj, Wi, Wj, gW0i, gW0j, W0i, W0j;
  Vector gradWi, gradWj, gradW0i, gradW0j;
  Vector deltagrad;

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (auto itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    const auto& nodeList = **itr;
    const auto  firstGhostNodei = nodeList.firstGhostNode();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  maxNumNeighbors = nodeList.maxNumNeighbors();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    // Get the work field for this NodeList.
    auto& workFieldi = nodeList.work();

    // Build the functor we use to compute the effective coupling between nodes.
    DamagedNodeCouplingWithFrags<Dimension> coupling(damage, gradDamage, H, fragIDs);

    // Check if we can identify a reference density.
    Scalar rho0 = 0.0;
    try {
      rho0 = dynamic_cast<const SolidEquationOfState<Dimension>&>(nodeList.equationOfState()).referenceDensity();
      // cerr << "Setting reference density to " << rho0 << endl;
    } catch(...) {
      // cerr << "BLAGO!" << endl;
    }

    // Iterate over the internal nodes in this NodeList.
    const auto ni = connectivityMap.numNodes(nodeListi);
#pragma omp parallel for                                \
  firstprivate(Ai, Aj,                                  \
               gradAi, gradAj, forceij, forceji,        \
               Bi, Bj,                                  \
               Ci, Cj,                                  \
               gradBi, gradBj,                          \
               gradCi, gradCj,                          \
               gWi, gWj, Wi, Wj, gW0i, gW0j, W0i, W0j,  \
               gradWi, gradWj, gradW0i, gradW0j,        \
               deltagrad)
    for (auto iItr = 0; iItr < ni; ++iItr) {
      const auto i = connectivityMap.ithNode(nodeListi, iItr);

      // Prepare to accumulate the time.
      const auto start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  si = entropy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto  mui = mu(nodeListi, i);
      Ai = A(nodeListi, i);
      gradAi = gradA(nodeListi, i);
      if (order != CRKOrder::ZerothOrder) {
        Bi = B(nodeListi, i);
        gradBi = gradB(nodeListi, i);
      }
      if (order == CRKOrder::QuadraticOrder) {
        Ci = C(nodeListi, i);
        gradCi = gradC(nodeListi, i);
      }
      const auto Hdeti = Hi.Determinant();
      const auto weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Ai > 0.0);
      CHECK(Hdeti > 0.0);
      CHECK(weighti > 0.0);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure(nodeListi, i);
      auto& viscousWorki = viscousWork(nodeListi, i);
      auto& pairAccelerationsi = pairAccelerations(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);
      auto& gradRhoi = gradRho(nodeListi, i);
      auto& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {

          // Loop over the neighbors.
#pragma vector always
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            ++ncalc;

            // Get the state for node j
            const auto& rj = position(nodeListj, j);
            const auto  mj = mass(nodeListj, j);
            const auto& vj = velocity(nodeListj, j);
            const auto  rhoj = massDensity(nodeListj, j);
            const auto  epsj = specificThermalEnergy(nodeListj, j);
            const auto  sj = entropy(nodeListj, j);
            const auto  Pj = pressure(nodeListj, j);
            const auto& Hj = H(nodeListj, j);
            const auto  cj = soundSpeed(nodeListj, j);
            Aj = A(nodeListj, j);
            gradAj = gradA(nodeListj, j);
            if (order != CRKOrder::ZerothOrder) {
              Bj = B(nodeListj, j);
              gradBj = gradB(nodeListj, j);
            }
            if (order == CRKOrder::QuadraticOrder) {
              Cj = C(nodeListj, j);
              gradCj = gradC(nodeListj, j);
            }
            const auto& Sj = S(nodeListj, j);
            const auto  Hdetj = Hj.Determinant();
            const auto  weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
            CHECK(mj > 0.0);
            CHECK(rhoj > 0.0);
            CHECK(Hdetj > 0.0);
            CHECK(weightj > 0.0);

            // Node displacement.
            const auto rij = ri - rj;
            const auto etai = Hi*rij;
            const auto etaj = Hj*rij;
            const auto etaMagi = etai.magnitude();
            const auto etaMagj = etaj.magnitude();
            CHECK(etaMagi >= 0.0);
            CHECK(etaMagj >= 0.0);
            const auto vij = vi - vj;

            // Symmetrized kernel weight and gradient.
            CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, CRKSPHHydroBase<Dimension>::correctionOrder(),  rij,  etai, Hi, Hdeti,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi);
            CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, CRKSPHHydroBase<Dimension>::correctionOrder(), -rij, -etaj, Hj, Hdetj, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj);
            deltagrad = gradWj - gradWi;
            const auto gradWSPHi = (Hi*etai.unitVector())*W.gradValue(etai.magnitude(), Hdeti);

            // Find the damaged pair weighting scaling.
            const auto fij = coupling(nodeListi, i, nodeListj, j);
            CHECK(fij >= 0.0 and fij <= 1.0);

            // Find the effective weights of i->j and j->i.
            // const auto wi = 2.0*weighti*weightj/(weighti + weightj);
            const auto wij = 0.5*(weighti + weightj);

            // Zero'th and second moment of the node distribution -- used for the
            // ideal H calculation.
            if (voidPoint(nodeListi, i) == 0 and voidPoint(nodeListj, j) == 0) {
              const auto fweightij = nodeListi == nodeListj ? 1.0 : mj*rhoi/(mi*rhoj);
              const auto rij2 = rij.magnitude2();
              const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
              weightedNeighborSumi += fweightij*std::abs(gWi);
              massSecondMomenti +=    fweightij*gradWSPHi.magnitude2()*thpt;
            }

            // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
            const auto QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                      ri, etai, vi, rhoi, ci, Hi,
                                      rj, etaj, vj, rhoj, cj, Hj);
            const auto Qaccij = (rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second).dot(deltagrad);
            // const auto workQij = 0.5*(vij.dot(Qaccij));
            const auto workQi = rhoj*rhoj*QPiij.second.dot(vij).dot(deltagrad);                // CRK
            const auto Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
            maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                                 // We need tighter timestep controls on the Q with CRK
            effViscousPressurei += wij * Qi * Wj;
            viscousWorki += 0.5*wij*wij/mi*workQi;

            // Velocity gradient.
            DvDxi -= wij*vij.dyad(gradWj);
            localDvDxi -= fij*wij*vij.dyad(gradWj);

            // Mass density gradient.
            gradRhoi += wij*(rhoj - rhoi)*gradWj;

            // We treat positive and negative pressures distinctly, so split 'em up.
            const auto Pposi = max(0.0, Pi),
                       Pnegi = min(0.0, Pi),
                       Pposj = max(0.0, Pj),
                       Pnegj = min(0.0, Pj);

            // Compute the stress tensors.
            SymTensor sigmai, sigmaj;
            if (nodeListi == nodeListj) {
              sigmai = Si - Pnegi*SymTensor::one;
              sigmaj = Sj - Pnegj*SymTensor::one;
            }

            // We decide between RK and CRK for the momentum and energy equations based on the surface condition.
            // Momentum
            forceij = (true ? // surfacePoint(nodeListi, i) <= 1 ? 
                       0.5*wij*wij*((Pposi + Pposj)*deltagrad - fij*(sigmai + sigmaj)*deltagrad + Qaccij) :                    // Type III CRK interpoint force.
                       mi*wij*(((Pposj - Pposi)*gradWj - fij*(sigmaj - sigmai)*gradWj)/rhoi + rhoi*QPiij.first.dot(gradWj)));  // RK
            DvDti -= forceij/mi;
            if (compatibleEnergy) pairAccelerationsi.push_back(-forceij/mi);

            // Energy
            DepsDti += (true ? // surfacePoint(nodeListi, i) <= 1 ?
                        0.5*wij*wij*(Pposj*vij.dot(deltagrad) + fij*sigmaj.dot(vij).dot(deltagrad) + workQi)/mi :               // CRK
                        wij*rhoi*QPiij.first.dot(vij).dot(gradWj));                                                             // RK, Q term only -- adiabatic portion added later

            // Estimate of delta v (for XSPH).
            XSPHDeltaVi -= fij*wij*Wj*vij;
          }
        }
      }
      const auto numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not this->compatibleEnergyEvolution() or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // For a surface point, add the RK thermal energy evolution.
      // if (surfacePoint(nodeListi, i) > 1) DepsDti += (Si - Pi*SymTensor::one).doubledot(DvDxi)/rhoi;

      // Get the time for pairwise interactions.
      const auto deltaTimePair = Timing::difference(start, Timing::currentTime())/max(size_t(1), ncalc);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = smoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                            ri,
                                                            DvDxi,
                                                            hmin,
                                                            hmax,
                                                            hminratio,
                                                            nPerh);
      Hideali = smoothingScaleMethod.newSmoothingScale(Hi,
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

      // Optionally use damage to ramp down stress on damaged material.
      const auto Di = (mDamageRelieveRubble ? 
                       max(0.0, min(1.0, damage(nodeListi, i).Trace() - 1.0)) :
                       0.0);
      // Hideali = (1.0 - Di)*Hideali + Di*mHfield0(nodeListi, i);
      // DHDti = (1.0 - Di)*DHDti + Di*(mHfield0(nodeListi, i) - Hi)*0.25/dt;

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - (deformation.Trace()/3.0)*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;

      // In the presence of damage, add a term to reduce the stress on this point.
      DSDti = (1.0 - Di)*DSDti - Di*Si*0.25/dt;

      // Time evolution of the mass density.
      DrhoDti = -rhoi*localDvDxi.Trace();

      // We also adjust the density evolution in the presence of damage.
      if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - 0.25/dt*Di*(rhoi - rho0);

      // If needed finish the total energy derivative.
      if (this->evolveTotalEnergy()) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Increment the work for i.
      worki += Timing::difference(start, Timing::currentTime());
    }
  }
}

}
