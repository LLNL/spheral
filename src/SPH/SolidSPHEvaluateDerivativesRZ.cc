namespace Spheral {
namespace SPHSpace {

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SolidSPHHydroBaseRZ::
evaluateDerivatives(const Dim<2>::Scalar time,
                    const Dim<2>::Scalar dt,
                    const DataBase<Dim<2> >& dataBase,
                    const State<Dim<2> >& state,
                    StateDerivatives<Dim<2> >& derivatives) const {

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();
  const auto& WG = this->GradKernel();
  const auto& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;
  const auto W0 = W(0.0, 1.0);
  const auto WQ0 = WQ(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto XSPH = this->XSPH();
  const auto damageRelieveRubble = this->damageRelieveRubble();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto STT = state.fields(SolidFieldNames::deviatoricStressTT, 0.0);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  const auto gradDamage = state.fields(SolidFieldNames::damageGradient, Vector::zero);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(STT.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(gradDamage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);

  const auto& Hfield0 = this->Hfield0();

  // Derivative FieldLists.
  auto rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto rhoSumCorrection = derivatives.fields(HydroFieldNames::massDensityCorrection, 0.0);
  auto viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto DSDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto DSTTDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStressTT, 0.0);
  CHECK(rhoSum.size() == numNodeLists);
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
  CHECK(rhoSumCorrection.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);
  CHECK(DSTTDt.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (compatibleEnergy) {
    auto nodeListi = 0;
    for (auto itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (auto i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Some scratch variables.
  Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
  Tensor QPiij, QPiji;
  SymTensor sigmai, sigmaj;

  // Start our big loop over all FluidNodeLists.
  auto nodeListi = 0;
  for (auto itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const auto& nodeList = **itr;
    const auto firstGhostNodei = nodeList.firstGhostNode();
    const auto hmin = nodeList.hmin();
    const auto hmax = nodeList.hmax();
    const auto hminratio = nodeList.hminratio();
    const auto maxNumNeighbors = nodeList.maxNumNeighbors();
    const auto nPerh = nodeList.nodesPerSmoothingScale();
    const auto Hmin = 1.0/hmin * SymTensor::one;

    // The scale for the tensile correction.
    const auto WnPerh = W(1.0/nPerh, 1.0);

    // Get the work field for this NodeList.
    auto& workFieldi = nodeList.work();

    // Build the functor we use to compute the effective coupling between nodes.
    DamagedNodeCouplingWithFrags<Dimension> coupling(damage, gradDamage, H, fragIDs);

    // Check if we can identify a reference density.
    auto rho0 = 0.0;
    try {
      rho0 = dynamic_cast<const SolidEquationOfState<Dimension>&>(nodeList.equationOfState()).referenceDensity();
      // cerr << "Setting reference density to " << rho0 << endl;
    } catch(...) {
      // cerr << "BLAGO!" << endl;
    }

    // Iterate over the internal nodes in this NodeList.
    for (auto iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const auto i = *iItr;

      // Prepare to accumulate the time.
      const auto start = Timing::currentTime();
      auto ncalc = 0;

      // Get the state for node i.
      const auto posi = position(nodeListi, i);
      const auto ri = abs(posi.y());
      const auto circi = 2.0*M_PI*ri;
      const auto mi = mass(nodeListi, i);
      const auto mRZi = mi/circi;
      const auto vi = velocity(nodeListi, i);
      const auto rhoi = massDensity(nodeListi, i);
      const auto epsi = specificThermalEnergy(nodeListi, i);
      const auto Pi = pressure(nodeListi, i);
      const auto Hi = H(nodeListi, i);
      const auto ci = soundSpeed(nodeListi, i);
      const auto omegai = omega(nodeListi, i);
      const auto Si = S(nodeListi, i);
      const auto STTi = STT(nodeListi, i);
      const auto mui = mu(nodeListi, i);
      const auto Hdeti = Hi.Determinant();
      const auto safeOmegai = safeInv(omegai, tiny);
      const auto fragIDi = fragIDs(nodeListi, i);
      const auto pTypei = pTypes(nodeListi, i);
      const auto zetai = abs((Hi*posi).y());
      const auto hri = ri*safeInv(zetai);
      const auto riInv = safeInv(ri, 0.25*hri);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure(nodeListi, i);
      auto& rhoSumCorrectioni = rhoSumCorrection(nodeListi, i);
      auto& viscousWorki = viscousWork(nodeListi, i);
      auto& pairAccelerationsi = pairAccelerations(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);
      auto& DSTTDti = DSTTDt(nodeListi, i);
      auto& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const auto firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const auto j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j.
              const auto posj = position(nodeListj, j);
              const auto rj = abs(posj.y());
              const auto circj = 2.0*M_PI*rj;
              const auto mj = mass(nodeListj, j);
              const auto mRZj = mj/circj;
              const auto vj = velocity(nodeListj, j);
              const auto rhoj = massDensity(nodeListj, j);
              const auto epsj = specificThermalEnergy(nodeListj, j);
              const auto Pj = pressure(nodeListj, j);
              const auto Hj = H(nodeListj, j);
              const auto cj = soundSpeed(nodeListj, j);
              const auto omegaj = omega(nodeListj, j);
              const auto Sj = S(nodeListj, j);
              const auto STTj = STT(nodeListj, j);
              const auto muj = mu(nodeListj, j);
              const auto Hdetj = Hj.Determinant();
              const auto safeOmegaj = safeInv(omegaj, tiny);
              const auto fragIDj = fragIDs(nodeListj, j);
              const auto pTypej = pTypes(nodeListj, j);
              const auto zetaj = abs((Hj*posj).y());
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Hdetj > 0.0);

              auto& rhoSumj = rhoSum(nodeListj, j);
              auto& DxDtj = DxDt(nodeListj, j);
              auto& DvDtj = DvDt(nodeListj, j);
              auto& DepsDtj = DepsDt(nodeListj, j);
              auto& DvDxj = DvDx(nodeListj, j);
              auto& localDvDxj = localDvDx(nodeListj, j);
              auto& Mj = M(nodeListj, j);
              auto& localMj = localM(nodeListj, j);
              auto& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              auto& effViscousPressurej = effViscousPressure(nodeListj, j);
              auto& rhoSumCorrectionj = rhoSumCorrection(nodeListj, j);
              auto& viscousWorkj = viscousWork(nodeListj, j);
              auto& pairAccelerationsj = pairAccelerations(nodeListj, j);
              auto& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
              auto& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              auto& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              auto& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Flag if this is a contiguous material pair or not.
              const auto sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

              // Flag if at least one particle is free (0).
              const auto freeParticle = (pTypei == 0 or pTypej == 0);

              // Node displacement.
              const auto xij = posi - posj;
              const auto etai = Hi*xij;
              const auto etaj = Hj*xij;
              const auto etaMagi = etai.magnitude();
              const auto etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);

              // Symmetrized kernel weight and gradient.
              std::tie(Wi, gWi) = W.kernelAndGradValue(etaMagi, Hdeti);
              std::tie(WQi, gWQi) = WQ.kernelAndGradValue(etaMagi, Hdeti);
              const auto Hetai = Hi*etai.unitVector();
              const auto gradWi = gWi*Hetai;
              const auto gradWQi = gWQi*Hetai;
              const auto gradWGi = WG.gradValue(etaMagi, Hdeti) * Hetai;

              std::tie(Wj, gWj) = W.kernelAndGradValue(etaMagj, Hdetj);
              std::tie(WQj, gWQj) = WQ.kernelAndGradValue(etaMagj, Hdetj);
              const auto Hetaj = Hj*etaj.unitVector();
              const auto gradWj = gWj*Hetaj;
              const auto gradWQj = gWQj*Hetaj;
              const auto gradWGj = WG.gradValue(etaMagj, Hdetj) * Hetaj;

              // Determine how we're applying damage.
              const auto fDeffij = coupling(nodeListi, i, nodeListj, j);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const auto fweightij = sameMatij ? 1.0 : mRZj*rhoi/(mRZi*rhoj);
              const auto xij2 = xij.magnitude2();
              const auto thpt = xij.selfdyad()*safeInvVar(xij2*xij2*xij2);
              weightedNeighborSumi +=     fweightij*abs(gWi);
              weightedNeighborSumj += 1.0/fweightij*abs(gWj);
              massSecondMomenti +=     fweightij*gradWi.magnitude2()*thpt;
              massSecondMomentj += 1.0/fweightij*gradWj.magnitude2()*thpt;

              // Contribution to the sum density (only if the same material).
              if (nodeListi == nodeListj) {
                rhoSumi += mRZj*Wi;
                rhoSumj += mRZi*Wj;
              }

              // Contribution to the sum density correction
              rhoSumCorrectioni += mRZj * WQi / rhoj ;
              rhoSumCorrectionj += mRZi * WQj / rhoi ;

              // Compute the pair-wise artificial viscosity.
              const auto vij = vi - vj;
              std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                              ri, etai, vi, rhoi, ci, Hi,
                                              rj, etaj, vj, rhoj, cj, Hj);
              const auto Qacci = 0.5*(QPiij*gradWQi);
              const auto Qaccj = 0.5*(QPiji*gradWQj);
              const auto workQi = vij.dot(Qacci);
              const auto workQj = vij.dot(Qaccj);
              const auto Qi = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
              const auto Qj = rhoj*rhoj*(QPiji.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);
              effViscousPressurei += mRZj*Qi*WQi/rhoj;
              effViscousPressurej += mRZi*Qj*WQj/rhoi;
              viscousWorki += mRZj*workQi;
              viscousWorkj += mRZi*workQj;

              // Damage scaling of negative pressures.
              const auto Peffi = (Pi > 0.0 ? Pi : fDeffij*Pi);
              const auto Peffj = (Pj > 0.0 ? Pj : fDeffij*Pj);

              // Compute the stress tensors.
              if (sameMatij) {
                sigmai = fDeffij*Si - Peffi*SymTensor::one;
                sigmaj = fDeffij*Sj - Peffj*SymTensor::one;
              } else {
                sigmai = -Peffi*SymTensor::one;
                sigmaj = -Peffj*SymTensor::one;
              }

              // Compute the tensile correction to add to the stress as described in 
              // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
              const auto fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
              const auto fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
              const auto Ri = fi*tensileStressCorrection(sigmai);
              const auto Rj = fj*tensileStressCorrection(sigmaj);
              sigmai += Ri;
              sigmaj += Rj;

              // Acceleration.
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              const auto sigmarhoi = safeOmegai*sigmai/(rhoi*rhoi);
              const auto sigmarhoj = safeOmegaj*sigmaj/(rhoj*rhoj);
              const auto deltaDvDt = sigmarhoi*gradWi + sigmarhoj*gradWj - Qacci - Qaccj;
              if (freeParticle) {
                DvDti += mRZj*deltaDvDt;
                DvDtj -= mRZi*deltaDvDt;
              }

              // Pair-wise portion of grad velocity.
              const auto deltaDvDxi = fDeffij*vij.dyad(gradWGi);
              const auto deltaDvDxj = fDeffij*vij.dyad(gradWGj);

              // Specific thermal energy evolution.
              DepsDti -= mRZj*(sigmarhoi.doubledot(deltaDvDxi.Symmetric()) - workQi);
              DepsDtj -= mRZi*(sigmarhoj.doubledot(deltaDvDxj.Symmetric()) - workQj);
              if (compatibleEnergy) {
                pairAccelerationsi.push_back( mRZj*deltaDvDt);
                pairAccelerationsj.push_back(-mRZi*deltaDvDt);
              }

              // Velocity gradient.
              DvDxi -= mRZj*deltaDvDxi;
              DvDxj -= mRZi*deltaDvDxj;
              if (sameMatij) {
                localDvDxi -= mRZj*deltaDvDxi;
                localDvDxj -= mRZi*deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (sameMatij or min(zetai, zetaj) < 1.0) {
                const auto wXSPHij = 0.5*(mRZi/rhoi*Wi + mRZj/rhoj*Wj);
                XSPHWeightSumi += wXSPHij;
                XSPHWeightSumj += wXSPHij;
                XSPHDeltaVi -= wXSPHij*vij;
                XSPHDeltaVj += wXSPHij*vij;
              }

              // Linear gradient correction term.
              Mi -= mRZj*xij.dyad(gradWGi);
              Mj -= mRZi*xij.dyad(gradWGj);
              if (sameMatij) {
                localMi -= mRZj*xij.dyad(gradWGi);
                localMj -= mRZi*xij.dyad(gradWGj);
              }
            }
          }
        }
      }
      const auto numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not this->compatibleEnergyEvolution() or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Get the time for pairwise interactions.
      const auto deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // Add the self-contribution to density sum.
      rhoSumi += mRZi*W0*Hdeti;
      rhoSumi /= circi;

      // Add the self-contribution to density sum correction.
      rhoSumCorrectioni += mRZi*WQ0*Hdeti/rhoi ;

      // Correct the effective viscous pressure.
      effViscousPressurei /= rhoSumCorrectioni ;

      // Finish the acceleration.
      const Vector deltaDvDti(Si(1,0)/rhoi*riInv,
                              (Si(1,1) - STTi)/rhoi*riInv);
      DvDti += deltaDvDti;
      pairAccelerationsi.push_back(deltaDvDti);

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      if (this->mCorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoi;
      }
      if (this->mCorrectVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      } else {
        localDvDxi /= rhoi;
      }

      // Evaluate the continuity equation.
      XSPHWeightSumi += Hdeti*mRZi/rhoi*W0;
      CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
      XSPHDeltaVi /= XSPHWeightSumi;
      const auto vri = vi.y(); // + XSPHDeltaVi.y();
      DrhoDti = -rhoi*(DvDxi.Trace() + vri*riInv);

      // Finish the specific thermal energy evolution.
      DepsDti += (STTi - Pi)/rhoi*vri*riInv;

      // If needed finish the total energy derivative.
      if (this->mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

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
                                                            posi,
                                                            DvDxi,
                                                            hmin,
                                                            hmax,
                                                            hminratio,
                                                            nPerh);
      Hideali = smoothingScaleMethod.newSmoothingScale(Hi,
                                                       posi,
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
      const auto Di = (damageRelieveRubble ? 
                       max(0.0, min(1.0, damage(nodeListi, i).Trace() - 1.0)) :
                       0.0);
      // Hideali = (1.0 - Di)*Hideali + Di*Hfield0(nodeListi, i);
      // DHDti = (1.0 - Di)*DHDti + Di*(Hfield0(nodeListi, i) - Hi)*0.25/dt;

      // We also adjust the density evolution in the presence of damage.
      if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - 0.25/dt*Di*(rhoi - rho0);

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto deformationTT = vi.y()*riInv;
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - ((deformation.Trace() + deformationTT)/3.0)*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;
      DSTTDti = 2.0*mui*(deformationTT - (deformation.Trace() + deformationTT)/3.0);

      // In the presence of damage, add a term to reduce the stress on this point.
      DSDti = (1.0 - Di)*DSDti - Di*Si*0.25/dt;
      DSTTDti = (1.0 - Di)*DSTTDti - Di*STTi*0.25/dt;

      // Increment the work for i.
      worki += Timing::difference(start, Timing::currentTime());

      // Now add the pairwise time for each neighbor we computed here.
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const auto firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
          auto& workFieldj = nodeLists[nodeListj]->work();
#pragma vector always
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const auto j = *jItr;
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

}
}
