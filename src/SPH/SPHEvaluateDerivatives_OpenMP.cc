namespace Spheral {

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  //static double totalLoopTime = 0.0;
  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;
  const double tiny = 1.0e-30;
  const Scalar W0 = W(0.0, 1.0);

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
  auto rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
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
  auto viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
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

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (auto itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const auto& nodeList = **itr;
    const auto firstGhostNodei = nodeList.firstGhostNode();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  maxNumNeighbors = nodeList.maxNumNeighbors();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    // The scale for the tensile correction.
    const auto  WnPerh = W(1.0/nPerh, 1.0);

    // Get the work field for this NodeList.
    auto& workFieldi = nodeList.work();

    const auto ni = connectivityMap.numNodes(nodeListi);
#pragma omp parallel for \
  private(Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj, QPiij, QPiji)
    for (auto iItr = 0; iItr < ni; ++iItr) {
      const auto i = connectivityMap.ithNode(nodeListi, iItr);

      // Prepare to accumulate the time.
      const auto start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& omegai = omega(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  safeOmegai = safeInv(omegai, tiny);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
      auto& normi = normalization(nodeListi, i);
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
      auto& viscousWorki = viscousWork(nodeListi, i);
      auto& pairAccelerationsi = pairAccelerations(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& worki = workFieldi(i);

      auto  maxvp = maxViscousPressurei;

      // Get the connectivity info for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {

// #ifdef _OPENMP

// #ifdef USE_UVM
//           #pragma omp target parallel for    \
//           reduction(max: maxvp) \
//           reduction(+: ncalc, weightedNeighborSumi, rhoSumi, normi,  \
//                   effViscousPressurei, viscousWorki, DepsDti, XSPHWeightSumi ) \
//           reduction(vecadd: DvDti, XSPHDeltaVi ) \
//           reduction(symtensadd: massSecondMomenti ) \
// 	  reduction(tensadd: Mi, localMi, DvDxi, localDvDxi) 
// // #else
// //                          // // rhoSumj, normj,          DepsDtj, effViscousPressurej, viscousWorkj, XSPHWeightSumj, weightedNeighborSumj) \
// //                               // DvDtj, XSPHDeltaVj )                                           \, massSecondMomentj
// //           #pragma omp parallel for                                                \
// //             default(shared)                                             \
// //             reduction(max: maxvp, maxViscousPressurei) \
// //             reduction(+: ncalc,  \
// //                          rhoSumi, normi, DrhoDti, DepsDti, effViscousPressurei, viscousWorki, XSPHWeightSumi, weightedNeighborSumi, worki) \
// //             reduction(vecadd: DvDti, XSPHDeltaVi) \
// //             reduction(symtensadd: massSecondMomenti )                \
// //             reduction(tensadd: Mi, localMi, DvDxi, localDvDxi)
// #endif

// #endif
// for (int jct=0; jct < nj; ++jct) {
//   const int j = *(jItr0+jct);

          // Loop over the neighbors.
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const auto j = *jItr;
            ncalc++; 

            // Get the state for node j
            const auto& rj = position(nodeListj, j);
            const auto& mj = mass(nodeListj, j);
            const auto& vj = velocity(nodeListj, j);
            const auto& rhoj = massDensity(nodeListj, j);
            const auto& epsj = specificThermalEnergy(nodeListj, j);
            const auto& Pj = pressure(nodeListj, j);
            const auto& Hj = H(nodeListj, j);
            const auto& cj = soundSpeed(nodeListj, j);
            const auto& omegaj = omega(nodeListj, j);
            const auto  Hdetj = Hj.Determinant();
            const auto  safeOmegaj = safeInv(omegaj, tiny);
            CHECK(mj > 0.0);
            CHECK(rhoj > 0.0);
            CHECK(Hdetj > 0.0);

            // Flag if this is a contiguous material pair or not.
            const bool sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

            // Node displacement.
            const auto rij = ri - rj;
            const auto etai = Hi*rij;
            const auto etaj = Hj*rij;
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

            std::tie(Wj, gWj) = W.kernelAndGradValue(etaMagj, Hdetj);
            std::tie(WQj, gWQj) = WQ.kernelAndGradValue(etaMagj, Hdetj);
            const auto Hetaj = Hj*etaj.unitVector();
            const auto gradWj = gWj*Hetaj;
            const auto gradWQj = gWQj*Hetaj;

            // Zero'th and second moment of the node distribution -- used for the
            // ideal H calculation.
            const auto fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
            const auto rij2 = rij.magnitude2();
            const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
            weightedNeighborSumi +=     fweightij*std::abs(gWi);
            massSecondMomenti +=     fweightij*gradWi.magnitude2()*thpt;

            // Contribution to the sum density.
            if (nodeListi == nodeListj) {
              rhoSumi += mj*Wi;
              normi += mi/rhoi*Wi;
            }

            // Compute the pair-wise artificial viscosity.
            const auto vij = vi - vj;
            std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                            ri, etai, vi, rhoi, ci, Hi,
                                            rj, etaj, vj, rhoj, cj, Hj);
            const auto Qacci = 0.5*(QPiij*gradWQi);
            const auto Qaccj = 0.5*(QPiji*gradWQj);
            const auto workQi = vij.dot(Qacci);
            const auto Qi = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
            maxViscousPressurei = max(maxViscousPressurei, Qi);
            effViscousPressurei += mj*Qi*WQi/rhoj;
            viscousWorki += mj*workQi;

            // Determine an effective pressure including a term to fight the tensile instability.
            //             const auto fij = epsTensile*pow(Wi/(Hdeti*WnPerh), nTensile);
            const auto fij = mEpsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
            const auto Ri = fij*(Pi < 0.0 ? -Pi : 0.0);
            const auto Rj = fij*(Pj < 0.0 ? -Pj : 0.0);
            const auto Peffi = Pi + Ri;
            const auto Peffj = Pj + Rj;

            // Acceleration.
            CHECK(rhoi > 0.0);
            CHECK(rhoj > 0.0);
            const auto Prhoi = safeOmegai*Peffi/(rhoi*rhoi);
            const auto Prhoj = safeOmegaj*Peffj/(rhoj*rhoj);
            const auto deltaDvDt = Prhoi*gradWi + Prhoj*gradWj + Qacci + Qaccj;
            DvDti += -mj*deltaDvDt;

            // Specific thermal energy evolution.
            // const Scalar workQij = 0.5*(mj*workQi + mi*workQj);
            DepsDti += mj*(Prhoi*vij.dot(gradWi) + workQi);
            if (mCompatibleEnergyEvolution) {
              pairAccelerationsi.push_back(-mj*deltaDvDt);
            }

            // Velocity gradient.
            const auto deltaDvDxi = mj*vij.dyad(gradWi);
            DvDxi += -deltaDvDxi; 
            if (sameMatij) {
              localDvDxi += -deltaDvDxi; 
            }

            // Estimate of delta v (for XSPH).
            if (mXSPH and (sameMatij)) {
              const auto wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
              XSPHWeightSumi += wXSPHij;
              XSPHDeltaVi -= wXSPHij*vij;
            }

            // Linear gradient correction term.
            Mi += -mj*rij.dyad(gradWi);
            if (sameMatij) {
              localMi += -mj*rij.dyad(gradWi);
            }
          }
	  //#ifdef _OPENMP
	  //totalLoopTime += omp_get_wtime() - time1;
          //printf("loop time %lf\n", totalLoopTime );
          //#endif
        }
      }
      const auto numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not mCompatibleEnergyEvolution or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Get the time for pairwise interactions.
      const auto deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);
      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;
      normi += mi/rhoi*W0*Hdeti;

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
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        XSPHWeightSumi += Hdeti*mi/rhoi*W0;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDt(nodeListi, i) = vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi);
      } else {
        DxDt(nodeListi, i) = vi;
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

}
