namespace Spheral {

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                          StateDerivatives<Dimension>& derivatives) const {


  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // Get the SlideSurfaces.
  auto& slides = this->slideSurface();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  const auto fullyDamagedThreshold=1.0e-2;
  const auto tiny = std::numeric_limits<double>::epsilon();
  const auto W0 = W(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto totalEnergy = this->evolveTotalEnergy();
  const auto epsDiffusionCoeff = this->specificThermalEnergyDiffusionCoefficient();
  const auto rhoStabilizeCoeff = this->densityStabilizationCoefficient();
  const auto surfaceForceCoeff = this->surfaceForceCoefficient();
  const auto xsphCoeff = this->xsphCoefficient();
  const auto XSPH = xsphCoeff > tiny;
  const auto diffuseEnergy = epsDiffusionCoeff>tiny and compatibleEnergy;
  const auto stabilizeDensity = rhoStabilizeCoeff>tiny;
  const auto alwaysAverageKernels = (mKernelAveragingMethod==KernelAveragingMethod::AlwaysAverageKernels);
  const auto averageInterfaceKernels = (mKernelAveragingMethod==KernelAveragingMethod::AverageInterfaceKernels);
  const auto constructHLLC = (mInterfaceMethod == InterfaceMethod::HLLCInterface);
  const auto activateConstruction = !(mInterfaceMethod == InterfaceMethod::NoInterface);
  const auto oneOverDimension = 1.0/Dimension::nDim;

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  numPairs = pairs.size();
  const auto  nPerh = nodeLists[0]->nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);

  // Get the state and derivative FieldLists.
  const auto interfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
  const auto interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  const auto interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto rawPressure = state.fields(FSIFieldNames::rawPressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  const auto yield = state.fields(SolidFieldNames::yieldStrength, 0.0);
  const auto invJ2 = state.fields(FSIFieldNames::inverseEquivalentDeviatoricStress, 0.0);

  CHECK(fragIDs.size()==numNodeLists);
  CHECK(interfaceFraction.size() == numNodeLists);
  CHECK(interfaceNormals.size() == numNodeLists);
  CHECK(interfaceSmoothness.size() == numNodeLists);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(rawPressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(K.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);
  CHECK(yieldStrength.size() == numNodeLists);
  CHECK(invJ2.size() == numNodeLists);

  // Derivative FieldLists.
  const auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  const auto  localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  const auto  DepsDx = derivatives.fields(FSIFieldNames::specificThermalEnergyGradient, Vector::zero);
  const auto  DPDx = derivatives.fields(FSIFieldNames::pressureGradient, Vector::zero);
  auto  smoothedInterfaceNormals = derivatives.fields(FSIFieldNames::smoothedInterfaceNormals,Vector::zero);
  auto  newInterfaceFraction = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + FSIFieldNames::interfaceFraction, 0.0);
  auto  newInterfaceSmoothness = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + FSIFieldNames::interfaceSmoothness, 0.0);
  auto  newInterfaceNormals = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + FSIFieldNames::interfaceNormals, Vector::zero);
  auto  normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto  DSDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto& pairDepsDt = derivatives.getAny(HydroFieldNames::pairWork, vector<Scalar>());
  
  CHECK(smoothedInterfaceNormals.size() == numNodeLists);
  CHECK(newInterfaceFraction.size() == numNodeLists);
  CHECK(newInterfaceSmoothness.size() == numNodeLists);
  CHECK(newInterfaceNormals.size() == numNodeLists);
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
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if(compatibleEnergy){
    pairAccelerations.resize(numPairs);
    pairDepsDt.resize(2*numPairs);
  }

  this->computeMCorrection(time,dt,dataBase,state,derivatives);

// Now we calculate  the hydro deriviatives
// Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, Wj, gWj, PLineari, PLinearj, epsLineari, epsLinearj;
    Tensor QPiij, QPiji;
    SymTensor sigmai, sigmaj;
    Vector sigmarhoi, sigmarhoj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto newInterfaceSmoothness_thread = newInterfaceSmoothness.threadCopy(threadStack);
    auto newInterfaceFraction_thread = newInterfaceFraction.threadCopy(threadStack);
    auto newInterfaceNormals_thread = newInterfaceNormals.threadCopy(threadStack);
    auto smoothedInterfaceNormals_thread = smoothedInterfaceNormals.threadCopy(threadStack);
    auto normalization_thread = normalization.threadCopy(threadStack);
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DrhoDt_thread = DrhoDt.threadCopy(threadStack);
    auto DSDt_thread = DSDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    
#pragma omp for
    for (auto kk = 0u; kk < numPairs; ++kk) {
      
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& interfaceNormalsi = interfaceNormals(nodeListi,i);
      const auto& interfaceSmoothnessi = interfaceSmoothness(nodeListi,i);
      const auto& interfaceFractioni = interfaceFraction(nodeListi,i);
      const auto& ri = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi,i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& rPi = rawPressure(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& pTypei = pTypes(nodeListi, i);
      const auto  fragIDi = fragIDs(nodeListi, i);
      const auto  Yi = yield(nodeListi, i);
      const auto  invJ2i = invJ2(nodeListi, i);
      const auto  voli = mi/rhoi;
      const auto  mui = max(mu(nodeListi,i),tiny);
      const auto  Ki = max(tiny,K(nodeListi,i))+4.0/3.0*mui;
      const auto  Hdeti = Hi.Determinant();
      
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      const auto& DepsDxi = DepsDx(nodeListi, i);
      const auto& DPDxi = DPDx(nodeListi, i);
      const auto& Mi = M(nodeListi, i);
      auto& normi = normalization_thread(nodeListi,i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& smoothedInterfaceNormalsi = smoothedInterfaceNormals_thread(nodeListi,i);
      auto& newInterfaceNormalsi = newInterfaceNormals_thread(nodeListi,i);
      auto& newInterfaceSmoothnessi = newInterfaceSmoothness_thread(nodeListi,i);
      auto& newInterfaceFractioni = newInterfaceFraction_thread(nodeListi,i);

      // Get the state for node j
      const auto& interfaceNormalsj = interfaceNormals(nodeListj,j);
      const auto& interfaceSmoothnessj = interfaceSmoothness(nodeListj,j);
      const auto& interfaceFractionj = interfaceFraction(nodeListj,j);
      const auto& rj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj,j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& rPj = rawPressure(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& Sj = S(nodeListj, j);
      const auto& pTypej = pTypes(nodeListj, j);
      const auto  fragIDj = fragIDs(nodeListj, j);
      const auto  Yj = yield(nodeListj, j);
      const auto  invJ2j = invJ2(nodeListj, j);
      const auto  volj = mj/rhoj;
      const auto  muj = max(mu(nodeListj,j),tiny);
      const auto  Kj = max(tiny,K(nodeListj,j))+4.0/3.0*muj;
      const auto  Hdetj = Hj.Determinant();

      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      const auto& DepsDxj = DepsDx(nodeListj, j);
      const auto& DPDxj = DPDx(nodeListj, j);
      const auto& Mj = M(nodeListj,j);
      auto& normj = normalization_thread(nodeListj,j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& smoothedInterfaceNormalsj = smoothedInterfaceNormals_thread(nodeListj,j);
      auto& newInterfaceNormalsj = newInterfaceNormals_thread(nodeListj,j);
      auto& newInterfaceSmoothnessj = newInterfaceSmoothness_thread(nodeListj,j);
      auto& newInterfaceFractionj = newInterfaceFraction_thread(nodeListj,j);

      // line of action
      const auto rij = ri - rj;
      const auto rhatij = rij.unitVector();
      
      // decoupling and boolean switches
      //-------------------------------------------------------
      // Flag if this is a contiguous material pair or not.
      const auto sameMatij =  (nodeListi == nodeListj and fragIDi==fragIDj);
      const auto differentMatij = !sameMatij; 
      const auto averageKernelij = ( (differentMatij and averageInterfaceKernels) or alwaysAverageKernels);

      // Flag if at least one particle is free (0).
      const auto freeParticle = (pTypei == 0 or pTypej == 0);
      
      // pairwise damage and nodal damage
      //const auto fDij = (sameMatij ? pairs[kk].f_couple : 0.0);
      const auto Di = max(0.0, min(1.0, damage(nodeListi, i).dot(rhatij).magnitude()));
      const auto Dj = max(0.0, min(1.0, damage(nodeListj, j).dot(rhatij).magnitude()));
      const auto fDi =  (sameMatij ? (1.0-Di)*(1.0-Di) : 0.0 );
      const auto fDj =  (sameMatij ? (1.0-Dj)*(1.0-Dj) : 0.0 );
      const auto fDij = (sameMatij ? pow(1.0-std::abs(Di-Dj),2.0) : 0.0 );

      // is Pmin being activated (Pmin->zero for material interfaces)
      const auto pLimiti = (sameMatij ? (Pi-rhoi*ci*ci*1e-5) : 0.0);
      const auto pLimitj = (sameMatij ? (Pj-rhoj*cj*cj*1e-5) : 0.0);
      const auto pminActivei = (rPi < pLimiti);
      const auto pminActivej = (rPj < pLimitj);
      
      // decoupling criteria 
      const auto isExpanding = (ri-rj).dot(vi-vj) > 0.0;
      const auto isFullyDamaged = (fDi<fullyDamagedThreshold) or (fDj<fullyDamagedThreshold);
      const auto isPastAdhesionThreshold = pminActivei or pminActivej;
      const auto decouple = isExpanding and (isFullyDamaged and isPastAdhesionThreshold);

      // do we need to construct our interface velocity?
      const auto constructInterface = (fDij < 1.0-fullyDamagedThreshold) and activateConstruction;
      const auto negligableShearWave = max(mui,muj) < 1.0e-5*min(Ki,Kj);

      // do we reduce our deviatoric stress
      const auto isTensile = (((Si+Sj)-(Pi+Pj)*SymTensor::one).dot(rhatij)).dot(rhatij) > 0;
      const auto damageReduceStress = isTensile or differentMatij;

      // Kernels
      //--------------------------------------
      const auto Hij = 0.5*(Hi+Hj);
      const auto etaij = Hij*rij;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagij = etaij.magnitude();
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagij >= 0.0);
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      const auto Hetai = Hi*etai.unitVector();
      const auto Hetaj = Hj*etaj.unitVector();
      auto gradWi = gWi*Hetai;
      auto gradWj = gWj*Hetaj;
      auto gradWiMi = gradWi;
      auto gradWjMj = gradWj;

      // average our kernels
      const auto gradWij = 0.5*(gradWi+gradWj);
      const auto Wij = 0.5*(Wi+Wj);
      if(averageKernelij){
        const auto gWij = 0.5*(gWi+gWj);
        const auto geoGradWij = (gradWi*gradWj)*safeInv(gradWij)
        const auto geoWij = gWi*gWj*safeInv(Wij)
        const auto geogWij = Wi*Wj*safeInv(gWij)
        Wi = geoWij;
        Wj = geoWij;
        gWi = geogWij;
        gWj = geoWij;
        gradWi = geoGradWij;
        gradWj = geoGradWij;
      }

      if(this->correctVelocityGradient()){
        gradWiMi = Mi.Transpose()*gradWi;
        gradWjMj = Mj.Transpose()*gradWj;
      }

      // interface normals 
      //-----------------------------------------------------------
      const auto fSij = ( sameMatij ? 1.0 : -1.0);         // direction parameter
      const auto Aij = fSij*voli*volj*(gradWiMi+gradWjMj); // "surface area vector"
      newInterfaceNormalsi -= Aij;
      newInterfaceNormalsj += Aij;
      if(interfaceFractioni > tiny and interfaceFractionj > tiny){
        smoothedInterfaceNormalsi += fSij*volj*interfaceFractionj*interfaceNormalsj*Wi;
        smoothedInterfaceNormalsj += fSij*voli*interfaceFractioni*interfaceNormalsi*Wj;
      }
      if (differentMatij){
        const auto alignment = max(fSij*interfaceNormalsi.dot(interfaceNormalsj),0.0);
        newInterfaceSmoothnessi += alignment*volj*Wij;
        newInterfaceSmoothnessj += alignment*voli*Wij;
        newInterfaceFractioni += volj*Wij;
        newInterfaceFractionj += voli*Wij;
      }

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.
      //---------------------------------------------------------------
      const auto rij2 = rij.magnitude2();
      const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
      weightedNeighborSumi += abs(gWi);
      weightedNeighborSumj += abs(gWj);
      massSecondMomenti += gradWi.magnitude2()*thpt;
      massSecondMomentj += gradWj.magnitude2()*thpt;

      if (!decouple){

        // Stress state
        //---------------------------------------------------------------
        const auto rhoij = 0.5*(rhoi+rhoj); 
        const auto cij = 0.5*(ci+cj); 
        const auto vij = vi - vj;

        // raw AV
        std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                        ri, etaij, vi, rhoij, cij, Hij,  
                                        rj, etaij, vj, rhoij, cij, Hij); 

        // slide correction
        if (slides.isSlideSurface(nodeListi,nodeListj)){
          const auto slideCorr = slides.slideCorrection(interfaceSmoothnessi,
                                                        interfaceSmoothnessj,
                                                        interfaceNormalsi,
                                                        interfaceNormalsj,
                                                        vi,
                                                        vj);
          QPiij *= slideCorr;
          QPiji *= slideCorr;
        }

        // save our max pressure from the AV for each node
        maxViscousPressurei = max(maxViscousPressurei, rhoi*rhoj * QPiij.diagonalElements().maxAbsElement());
        maxViscousPressurej = max(maxViscousPressurej, rhoi*rhoj * QPiji.diagonalElements().maxAbsElement());

        // stress tensor
//        {
          // apply yield pairwise 
          const auto Yij = std::max(0.0,std::min(Yi,Yj));
          const auto fYieldi = std::min(Yij*invJ2i,1.0);
          const auto fYieldj = std::min(Yij*invJ2j,1.0);
          const auto Seffi = (damageReduceStress ? fDij : 1.0) * fYieldi * Si;
          const auto Seffj = (damageReduceStress ? fDij : 1.0) * fYieldj * Sj;
          const auto Peffi = (differentMatij ? max(Pi,0.0) : Pi);
          const auto Peffj = (differentMatij ? max(Pj,0.0) : Pj);
          sigmai = Seffi - Peffi * SymTensor::one;
          sigmaj = Seffj - Peffj * SymTensor::one;
//        }

        // Compute the tensile correction to add to the stress as described in 
        // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
        const auto fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
        const auto fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
        const auto Ri = fi*tensileStressCorrection(sigmai);
        const auto Rj = fj*tensileStressCorrection(sigmaj);
        sigmai += Ri;
        sigmaj += Rj;

        // accelerations
        //---------------------------------------------------------------
        const auto rhoirhoj = 1.0/(rhoi*rhoj);
        const auto sf = (sameMatij ? 1.0 : 1.0 + surfaceForceCoeff*abs((rhoi-rhoj)/(rhoi+rhoj+tiny)));
        sigmarhoi = sf*(rhoirhoj*sigmai-0.5*QPiij)*gradWiMi;
        sigmarhoj = sf*(rhoirhoj*sigmaj-0.5*QPiji)*gradWjMj;

        if (averageKernelij){
          const auto sigmarhoij = 0.5*(sigmarhoi+sigmarhoj);
          sigmarhoi = sigmarhoij;
          sigmarhoj = sigmarhoij;
        }
      
        const auto deltaDvDt = sigmarhoi + sigmarhoj;

        if (freeParticle) {
          DvDti += mj*deltaDvDt;
          DvDtj -= mi*deltaDvDt;
        } 
      
        // Velocity Gradient
        //-----------------------------------------------------------
        linearReconstruction(ri,rj,Pi,Pj,DPDxi,DPDxj,PLineari,PLinearj);
        // construct our interface velocity 
        auto vstar = 0.5*(vi+vj);
        
        if (constructInterface){
          
          // components
          const auto ui = vi.dot(rhatij);
          const auto uj = vj.dot(rhatij);
          const auto wi = vi - ui*rhatij;
          const auto wj = vj - uj*rhatij;
          
          // weights weights
          const auto Ci =  (constructHLLC ? std::sqrt(rhoi*Ki)  : Ki  ) + tiny;
          const auto Cj =  (constructHLLC ? std::sqrt(rhoj*Kj)  : Kj  ) + tiny;
          const auto Csi = (constructHLLC ? std::sqrt(rhoi*mui) : mui ) + tiny;
          const auto Csj = (constructHLLC ? std::sqrt(rhoj*muj) : muj ) + tiny;
          const auto CiCjInv = safeInv(Ci+Cj,tiny);
          const auto CsiCsjInv = safeInv(Csi+Csj,tiny);

          // weights
          const auto weightUi = max(0.0, min(1.0, Ci*CiCjInv));
          const auto weightUj = 1.0 - weightUi;
          const auto weightWi = (negligableShearWave ? weightUi : max(0.0, min(1.0, Csi*CsiCsjInv )) );
          const auto weightWj = 1.0 - weightWi;

          // interface velocity
          const auto ustar = weightUi*ui + weightUj*uj + (constructHLLC ? (PLinearj - PLineari)*CiCjInv : 0.0); 
          const auto wstar = weightWi*wi + weightWj*wj;// + (constructHLLC ? (Seffi - Seffj).dot(rhatij)*CsiCsjInv : Vector::zero); ;
          vstar = fDij * vstar + (1.0-fDij)*(ustar*rhatij + wstar);

          

          // const auto weightUi = max(0.0, min(1.0, Ci/(Ci+Cj)));
          // const auto weightUj = 1.0 - weightUi;
          // const auto weightWi = (negligableShearWave ? weightUi : max(0.0, min(1.0, Csi/(Csi+Csj) )) );
          // const auto weightWj = 1.0 - weightWi;

          // // get our eff pressure
          // const auto ustar = weightUi*ui + weightUj*uj; 
          // const auto wstar = weightWi*wi + weightWj*wj + fDij*(Si - Sj).dot(rhatij)/;
          // vstar = fDij * vstar + (1.0-fDij)*(ustar*rhatij + wstar);
  
        }

        // local velocity gradient for DSDt
        if (sameMatij) {
          localDvDxi -=  2.0*volj*((vi-vstar).dyad(gradWi));
          localDvDxj -=  2.0*voli*((vstar-vj).dyad(gradWj)); 
        }
        
        // diffuse to stabilize things
        if (stabilizeDensity and (ci>tiny and cj>tiny)){
          
          //const auto cFactor = 1.0 + max(min( (vi-vj).dot(rhatij)/max(cij,tiny), 0.0), -1.0);
          //const auto effCoeff = (differentMatij ? 1.0 : rhoStabilizeCoeff*cFactor);
          vstar +=  (1-fDij) * rhoStabilizeCoeff * rhatij * cij * min(max((PLinearj-PLineari)/(Ki + Kj),-0.25),0.25);
        }

        // global velocity gradient
        DvDxi -= 2.0*volj*(vi-vstar).dyad(gradWiMi);
        DvDxj -= 2.0*voli*(vstar-vj).dyad(gradWjMj);

        // energy conservation
        // ----------------------------------------------------------
        const auto deltaDepsDti = 2.0*sigmarhoi.dot(vi-vstar);
        const auto deltaDepsDtj = 2.0*sigmarhoj.dot(vstar-vj);

        DepsDti -= mj*deltaDepsDti;
        DepsDtj -= mi*deltaDepsDtj;

        if(compatibleEnergy){
          pairAccelerations[kk] = - deltaDvDt;
          pairDepsDt[2*kk]   = - deltaDepsDti; 
          pairDepsDt[2*kk+1] = - deltaDepsDtj;
        }
        
        // thermal diffusion
        //-----------------------------------------------------------
        if (sameMatij and diffuseEnergy){
          linearReconstruction(ri,rj,epsi,epsj,DepsDxi,DepsDxj,epsLineari,epsLinearj);
          const auto cijEff = max(min(cij + (vi-vj).dot(rhatij), cij),0.0);
          const auto diffusion =  epsDiffusionCoeff*cijEff*(epsLineari-epsLinearj)*etaij.dot(gradWij)/(rhoij*etaMagij*etaMagij+tiny);
          pairDepsDt[2*kk]   += diffusion; 
          pairDepsDt[2*kk+1] -= diffusion;
        }

        // normalization 
        //-----------------------------------------------------------
        normi += volj*Wi;
        normj += voli*Wj;

        // XSPH -- we use this to handle tensile instability here
        //-----------------------------------------------------------
        if (sameMatij and XSPH) {
          const auto fxsph  = (min(Pi,Pj) < 0.0 ? 1.0 : 0.0);
          XSPHWeightSumi += fxsph*volj*Wi;
          XSPHWeightSumj += fxsph*voli*Wj;
          XSPHDeltaVi -= volj*Wi*(vi-vstar);
          XSPHDeltaVj -= voli*Wj*(vj-vstar);
        }

      } // if damageDecouple 
    } // loop over pairs
    threadReduceFieldLists<Dimension>(threadStack);
  } // OpenMP parallel region


  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& mui = mu(nodeListi, i);
      const auto& interfaceFractioni = interfaceFraction(nodeListi,i);
      const auto& interfaceNormalsi = interfaceNormals(nodeListi,i);
      const auto  Hdeti = Hi.Determinant();
      
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      const auto& DvDti = DvDt(nodeListi,i);
      const auto& localMi = localM(nodeListi, i);
      auto& normi = normalization(nodeListi,i);
      auto& DepsDti = DepsDt(nodeListi,i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);
      auto& newInterfaceSmoothnessi = newInterfaceSmoothness(nodeListi,i);
      auto& newInterfaceFractioni = newInterfaceFraction(nodeListi,i);
      auto& newInterfaceNormalsi = newInterfaceNormals(nodeListi,i);
      auto& smoothedInterfaceNormalsi = smoothedInterfaceNormals(nodeListi,i);

      // finish our normalization
      normi += Hdeti*mi/rhoi*W0;

      // finish our interface fields.
      newInterfaceSmoothnessi = min(1.0,max(0.0,newInterfaceSmoothnessi/max(newInterfaceFractioni,tiny)));
      smoothedInterfaceNormalsi +=  interfaceFractioni * Hdeti*mi/rhoi*W0 * interfaceNormalsi;
      if (newInterfaceFractioni > tiny){
        const auto normalSmoothFraction = min(0.9,max(0.1, 20.0*min(interfaceFractioni,0.05)));
        newInterfaceNormalsi = (normalSmoothFraction * newInterfaceNormalsi.unitVector()+
                                (1.0-normalSmoothFraction) * smoothedInterfaceNormalsi.unitVector()).unitVector();
      }else{
        newInterfaceNormalsi = Vector::zero;
      }

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;
 
    
      DrhoDti -=  rhoi*DvDxi.Trace();

      if (totalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      DxDti = vi;
      if (XSPH) {
        XSPHWeightSumi += Hdeti*mi/rhoi*W0;
        XSPHWeightSumi /= normi;
        //CHECK(XSPHWeightSumi >= 0.0);
        DxDti += xsphCoeff*XSPHDeltaVi/max(normi,tiny);
      }

    
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

      if(this->correctVelocityGradient()) localDvDxi = localDvDxi*localMi;

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - (deformation.Trace()*oneOverDimension)*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti += spinCorrection + 2.0*mui*deviatoricDeformation;
      
    } //loop-nodes
  } //loop-nodeLists
} // evaluateDerivatives method


//------------------------------------------------------------------------------
// EvalDerivs subroutine for spatial derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
computeMCorrection(const typename Dimension::Scalar /*time*/,
                   const typename Dimension::Scalar /*dt*/,
                   const DataBase<Dimension>& dataBase,
                   const State<Dimension>& state,
                         StateDerivatives<Dimension>& derivatives) const {

  // The kernels and such.
  const auto& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  const auto alwaysAverageKernels = (mKernelAveragingMethod==KernelAveragingMethod::AlwaysAverageKernels);
  const auto averageInterfaceKernels = (mKernelAveragingMethod==KernelAveragingMethod::AverageInterfaceKernels);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  numPairs = pairs.size();

  // Get the state and derivative FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DepsDx = derivatives.fields(FSIFieldNames::specificThermalEnergyGradient, Vector::zero);
  auto  DPDx = derivatives.fields(FSIFieldNames::pressureGradient, Vector::zero);
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  
  CHECK(DepsDx.size() == numNodeLists);
  CHECK(DPDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);

#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    auto DPDx_thread = DPDx.threadCopy(threadStack);
    auto DepsDx_thread = DepsDx.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < numPairs; ++kk) {

      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& fragIDi = fragIDs(nodeListi,i);
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      // Get the state for node j
      const auto& fragIDj = fragIDs(nodeListj,j);
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& DPDxi = DPDx_thread(nodeListi, i);
      auto& DPDxj = DPDx_thread(nodeListj, j);
      auto& DepsDxi = DepsDx_thread(nodeListi, i);
      auto& DepsDxj = DepsDx_thread(nodeListj, j);
      auto& localMi = localM_thread(nodeListi,i);
      auto& localMj = localM_thread(nodeListj,j);
      auto& Mi = M_thread(nodeListi,i);
      auto& Mj = M_thread(nodeListj,j);

      const auto rij = ri - rj;
      const auto Pij = Pi - Pj;
      const auto epsij = epsi - epsj;

      // logic
      //---------------------------------------
      const auto sameMatij = (nodeListi == nodeListj and fragIDi == fragIDj);
      const auto differentMatij = (nodeListi!=nodeListj);
      const auto averageKernelij = ( (differentMatij and averageInterfaceKernels) or alwaysAverageKernels);

      // Kernels
      //--------------------------------------
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      const auto gWi = W.gradValue(etaMagi, Hdeti);
      const auto gWj = W.gradValue(etaMagj, Hdetj);
      
      const auto Hetai = Hi*etai.unitVector();
      const auto Hetaj = Hj*etaj.unitVector();
      
      auto gradWi = gWi*Hetai;
      auto gradWj = gWj*Hetaj;
      
      //Wi & Wj --> Wij for interface better agreement DrhoDt and DepsDt
      if(averageKernelij){
        const auto gradWij = 0.5*(gradWi+gradWj);
        const auto geoGradWij = (gradWi*gradWj)*safeInv(gradWij)
        gradWi = geoGradWij;
        gradWj = geoGradWij;
      }

      gradWi *= mj/rhoj;
      gradWj *= mi/rhoi;

      // spatial gradients and correction
      //---------------------------------------------------------------
      const auto deltaRi = rij.dyad(gradWi);
      const auto deltaRj = rij.dyad(gradWj);

      Mi -= deltaRi;
      Mj -= deltaRj;

      DPDxi -= Pij*gradWi;
      DPDxj -= Pij*gradWj;

      if(sameMatij){
        localMi -=  deltaRi;
        localMj -=  deltaRj;
        DepsDxi -= epsij*gradWi;
        DepsDxj -= epsij*gradWj;
      }
    } // loop over pairs
      // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
  }   // OpenMP parallel region

   
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto& nodeList = mass[nodeListi]->nodeList();
      const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
        const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
        auto& Mi = M(nodeListi, i);
        auto& localMi = localM(nodeListi, i);
        auto& DepsDxi = DepsDx(nodeListi, i);
        auto& DPDxi = DPDx(nodeListi, i);

        const auto Mdeti = Mi.Determinant();
        const auto goodM = ( Mdeti > 1.0e-2 and numNeighborsi > Dimension::pownu(2));
        Mi =  (goodM ? Mi.Inverse() : Tensor::one);
        
        const auto localMdeti = localMi.Determinant();
        const auto goodLocalM = ( localMdeti > 1.0e-2 and numNeighborsi > Dimension::pownu(2));
        localMi =  (goodLocalM ? localMi.Inverse() : Tensor::one);

        DPDxi = Mi.Transpose()*DPDxi;
        DepsDxi = localMi.Transpose()*DepsDxi;

      } // for each node
    }   // for each nodelist

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
      (*boundaryItr)->applyFieldListGhostBoundary(M);
      (*boundaryItr)->applyFieldListGhostBoundary(DPDx);
      (*boundaryItr)->applyFieldListGhostBoundary(DepsDx);
    }
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

} // method 

} // spheral namespace
