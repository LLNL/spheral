namespace Spheral {

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
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
  const auto tiny = std::numeric_limits<double>::epsilon();
  const auto W0 = W(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto totalEnergy = this->evolveTotalEnergy();
  const auto damageRelieveRubble = this->damageRelieveRubble();
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
  //const auto oneOverDimension = 1.0/Dimension::nDim;

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero);
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
  CHECK(S.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DepsDx = derivatives.fields(FSIFieldNames::specificThermalEnergyGradient, Vector::zero);
  auto  DPDx = derivatives.fields(FSIFieldNames::pressureGradient, Vector::zero);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
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

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Size up the pair-wise accelerations before we start.
  if(compatibleEnergy){
    pairAccelerations.resize(npairs);
    pairDepsDt.resize(2*npairs);
  }

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);

//M corr needs to be calculated beforehand 
//to be consistently applied to the acceleration
//and the time derivative of internal energy
if(this->correctVelocityGradient()){


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
    for (auto kk = 0u; kk < npairs; ++kk) {

      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
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

      // logic
      //---------------------------------------
      const auto sameMatij = (nodeListi == nodeListj);
      const auto differentMatij = (nodeListi!=nodeListj);
      const auto averageKernelij = ( (differentMatij and averageInterfaceKernels) or alwaysAverageKernels);
      
      // Kernels
      //--------------------------------------
      const auto rij = ri - rj;
      const auto Pij = Pi - Pj;
      const auto epsij = epsi - epsj;

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
        gradWi = gradWij;
        gradWj = gradWij;
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

      } 
    }

  //localM.Zero();

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

} // if correctVelocityGradient



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
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DrhoDt_thread = DrhoDt.threadCopy(threadStack);
    auto DSDt_thread = DSDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    //auto localM_thread = localM.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    
#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi,i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& pTypei = pTypes(nodeListi, i);
      const auto  voli = mi/rhoi;
      const auto  mui = max(mu(nodeListi,i),tiny);
      const auto  Ki = max(tiny,rhoi*ci*ci)+4.0/3.0*mui;
      auto  Hdeti = Hi.Determinant();
      //const auto fragIDi = fragIDs(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      const auto& DepsDxi = DepsDx(nodeListi, i);
      const auto& DPDxi = DPDx(nodeListi, i);
      const auto& Mi = M(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      //auto& localMi = localM_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      
      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj,j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& Sj = S(nodeListj, j);
      const auto& pTypej = pTypes(nodeListj, j);
      const auto  volj = mj/rhoj;
      const auto  muj = max(mu(nodeListj,j),tiny);
      const auto  Kj = max(tiny,rhoj*cj*cj)+4.0/3.0*muj;
      auto  Hdetj = Hj.Determinant();
      //const auto fragIDj = fragIDs(nodeListj, j);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      const auto& DepsDxj = DepsDx(nodeListj, j);
      const auto& DPDxj = DPDx(nodeListj, j);
      const auto& Mj = M(nodeListj,j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      //auto& localMj = localM_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      
      // Flag if this is a contiguous material pair or not.
      const auto sameMatij =  (nodeListi == nodeListj);// and fragIDi == fragIDj); 
      const auto differentMatij = !sameMatij;
      const auto averageKernelij = ( (differentMatij and averageInterfaceKernels) or alwaysAverageKernels);

      // Flag if at least one particle is free (0).
      const auto freeParticle = (pTypei == 0 or pTypej == 0);

      // we'll need a couple damage defs
      const auto rij = ri - rj;
      const auto rhatij = rij.unitVector();
      //const auto fDij = (sameMatij ? pairs[kk].f_couple : 0.0);
      const auto  Di = max(0.0, min(1.0, damage(nodeListi, i).dot(rhatij).magnitude()));
      const auto  Dj = max(0.0, min(1.0, damage(nodeListj, j).dot(rhatij).magnitude()));
      const auto fSij = (sameMatij ? 1.0-abs(Di-Dj) : 0.0);
      const auto fDi =  (sameMatij ? (1.0-Di)*(1.0-Di) : 0.0 );
      const auto fDj =  (sameMatij ? (1.0-Dj)*(1.0-Dj) : 0.0 );

      // Decoupling
      //-------------------------------------------------------
      // we need to test if these nodes are allowed to interact
      const auto isExpanding = (ri-rj).dot(vi-vj) > 0.0;
      const auto cantSupportTension = (fDi<0.01) or (fDj<0.01);
      const auto isInTension = (Pi<0.0) or (Pj<0.0);

      const auto decouple =  isExpanding and (cantSupportTension and isInTension);

      const auto constructInterface = (fSij < 0.99) and activateConstruction;
      const auto negligableShearWave = max(mui,muj) < 1.0e-5*min(Ki,Kj);

      if (!decouple){

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
        std::tie(Wi, gWi) = W.kernelAndGradValue(etaMagi, Hdeti);
        std::tie(Wj, gWj) = W.kernelAndGradValue(etaMagj, Hdetj);
        const auto Hetai = Hi*etai.unitVector();
        const auto Hetaj = Hj*etaj.unitVector();
        auto gradWi = gWi*Hetai;
        auto gradWj = gWj*Hetaj;
        auto gradWiMi = gradWi;
        auto gradWjMj = gradWj;

        // average our kernels
        const auto gradWij = 0.5*(gradWi+gradWj);
        if(averageKernelij){
          const auto Wij = 0.5*(Wi+Wj);
          const auto gWij = 0.5*(gWi+gWj);
          Wi = Wij;
          Wj = Wij;
          gWi = gWij;
          gWj = gWij;
          gradWi = gradWij;
          gradWj = gradWij;
        }

        if(this->correctVelocityGradient()){
         gradWiMi = Mi.Transpose()*gradWi;
         gradWjMj = Mj.Transpose()*gradWj;
        }

        // Zero'th and second moment of the node distribution -- used for the
        // ideal H calculation.
        //---------------------------------------------------------------
        const auto rij2 = rij.magnitude2();
        const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
        weightedNeighborSumi += abs(gWi);
        weightedNeighborSumj += abs(gWj);
        massSecondMomenti += gradWiMi.magnitude2()*thpt;
        massSecondMomentj += gradWjMj.magnitude2()*thpt;

        // Stress state
        //---------------------------------------------------------------
        const auto rhoij = 0.5*(rhoi+rhoj); 
        const auto cij = 0.5*(ci+cj); 
        const auto vij = vi - vj;
        
        // artificial viscosity
        std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                        ri, etaij, vi, rhoij, cij, Hij,  
                                        rj, etaij, vj, rhoij, cij, Hij); 
        
        const auto slideCorrection = slides.slideCorrection(nodeListi, i, nodeListj, j,vi,vj);

        QPiij *= slideCorrection;
        QPiji *= slideCorrection;

        maxViscousPressurei = max(maxViscousPressurei, rhoi*rhoj * QPiij.diagonalElements().maxAbsElement());
        maxViscousPressurej = max(maxViscousPressurej, rhoi*rhoj * QPiji.diagonalElements().maxAbsElement());

        // stress tensor
        const auto Peffi = (Pi<0.0 ? fDi : 1.0) * Pi;
        const auto Peffj = (Pj<0.0 ? fDj : 1.0) * Pj;

        sigmai = fDi * Si - Peffi * SymTensor::one;
        sigmaj = fDj * Sj - Peffj * SymTensor::one;

        // Compute the tensile correction to add to the stress as described in 
        // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
        const auto fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
        const auto fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
        const auto Ri = fi*tensileStressCorrection(sigmai);
        const auto Rj = fj*tensileStressCorrection(sigmaj);
        sigmai += Ri;
        sigmaj += Rj;

        if (averageKernelij){
          const auto sigmaij = 0.5*(sigmai+sigmaj);
          sigmai = sigmaij;
          sigmaj = sigmaij;
        }

        // accelerations
        //---------------------------------------------------------------
        const auto rhoirhoj = 1.0/(rhoi*rhoj);
        const auto sf = (sameMatij ? 1.0 : 1.0 + surfaceForceCoeff*abs((rhoi-rhoj)/(rhoi+rhoj+tiny)));
        sigmarhoi = sf*((rhoirhoj*sigmai-0.5*QPiij))*gradWiMi;
        sigmarhoj = sf*((rhoirhoj*sigmaj-0.5*QPiji))*gradWjMj;
      
        const auto deltaDvDt = sigmarhoi + sigmarhoj;

        if (freeParticle) {
          DvDti += mj*deltaDvDt;
          DvDtj -= mi*deltaDvDt;
        } 
      
        // Velocity Gradient
        //-----------------------------------------------------------

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

          const auto weightUi = max(0.0, min(1.0, Ci/(Ci+Cj)));
          const auto weightUj = 1.0 - weightUi;
          const auto weightWi = (negligableShearWave ? weightUi : max(0.0, min(1.0, Csi/(Csi+Csj) )) );
          const auto weightWj = 1.0 - weightWi;

          // get our eff pressure
          const auto ustar = weightUi*ui + weightUj*uj; 
          const auto wstar = weightWi*wi + weightWj*wj;
          vstar = fSij * vstar + (1.0-fSij)*(ustar*rhatij + wstar);
  
        }

        // local velocity gradient for DSDt
        if (sameMatij) {
          localDvDxi -=  2.0*volj*((vi-vstar).dyad(gradWi));
          localDvDxj -=  2.0*voli*((vstar-vj).dyad(gradWj)); 
        }

        // diffuse to stabilize things
        if (stabilizeDensity and (ci>tiny and cj>tiny)){
          linearReconstruction(ri,rj,Pi,Pj,DPDxi,DPDxj,PLineari,PLinearj);
          const auto cFactor = 1.0 + max(min( (vi-vj).dot(rhatij)/max(cij,tiny), 0.0), -1.0);
          const auto effCoeff = (differentMatij ? 1.0 : rhoStabilizeCoeff*cFactor);
          vstar += effCoeff * rhatij * cij * min(max((PLinearj-PLineari)/(Ki + Kj),-0.25),0.25);
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
        
        // diffusion
        //-----------------------------------------------------------
        if (sameMatij and diffuseEnergy){
          linearReconstruction(ri,rj,epsi,epsj,DepsDxi,DepsDxj,epsLineari,epsLinearj);
          const auto cijEff = max(min(cij + (vi-vj).dot(rhatij), cij),0.0);
          const auto diffusion =  epsDiffusionCoeff*cijEff*(epsLineari-epsLinearj)*etaij.dot(gradWij)/(rhoij*etaMagij*etaMagij+tiny);
          pairDepsDt[2*kk]   += diffusion; 
          pairDepsDt[2*kk+1] -= diffusion;
        }


        // XSPH
        //-----------------------------------------------------------
        if (XSPH) {
          XSPHWeightSumi += volj*Wi;
          XSPHWeightSumj += voli*Wj;
          XSPHDeltaVi -= 2.0*volj*Wi*(vi-vstar);
          XSPHDeltaVj -= 2.0*voli*Wj*(vj-vstar);
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
      const auto  Hdeti = Hi.Determinant();
      //const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      const auto& DvDti = DvDt(nodeListi,i);
      const auto& localMi = localM(nodeListi, i);
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
      
      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;
 
      if (totalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      DrhoDti -=  rhoi*DvDxi.Trace();

      DxDti = vi;
      if (XSPH) {
        CHECK(XSPHWeightSumi >= 0.0);
        XSPHWeightSumi += Hdeti*mi/rhoi*W0 + tiny;
        DxDti += xsphCoeff*XSPHDeltaVi/XSPHWeightSumi;
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

      // kill S when fully damaged
      const auto Di = (damageRelieveRubble ? 
                       max(0.0, min(1.0, damage(nodeListi, i).Trace()/Dimension::nDim)) :
                       0.0);

      //const auto localMdeti = localMi.Determinant();
      //const auto goodLocalM = ( localMdeti > 1.0e-2 and numNeighborsi > Dimension::pownu(2));
      //localMi =  (goodLocalM ? localMi.Inverse() : Tensor::one);
      localDvDxi = localDvDxi*localMi;

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - (deformation.Trace()/Dimension::nDim)*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti += spinCorrection + 2.0*mui*deviatoricDeformation;

      // In the presence of damage, add a term to reduce the stress on this point.
      if(Di>0.99) DSDti = (1.0 - Di)*DSDti - 0.125/dt*Di*Si;
      
    } //loop-nodes
  } //loop-nodeLists
} // evaluateDerivatives method


}
