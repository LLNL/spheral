namespace Spheral {
//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                          StateDerivatives<Dimension>& derivatives) const {

  const auto& riemannSolver = this->riemannSolver();

  const auto& smoothingScale = this->smoothingScaleMethod();
  
  // A few useful constants we'll use in the following loop.
  const auto tiny = std::numeric_limits<Scalar>::epsilon();
  const auto xsph = this->XSPH();
  const auto epsTensile = this->epsilonTensile();
  //const auto epsDiffusionCoeff = this->specificThermalEnergyDiffusionCoefficient();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  //const auto totalEnergy = this->evolveTotalEnergy();
  const auto gradType = this->gradientType();
  //const auto correctVelocityGradient = this->correctVelocityGradient();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  const auto  numNodeLists = nodeLists.size();
  const auto  nPerh = nodeLists[0]->nodesPerSmoothingScale();

  // kernel
  const auto& W = this->kernel();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  const auto  W0 = W(0.0, 1.0);

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto nodalVelocity = state.fields(GSPHFieldNames::nodalVelocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto riemannDpDx = state.fields(GSPHFieldNames::RiemannPressureGradient,Vector::zero);
  const auto riemannDvDx = state.fields(GSPHFieldNames::RiemannVelocityGradient,Tensor::zero);
  
  CHECK(nodalVelocity.size() == numNodeLists);  
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(volume.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(riemannDpDx.size() == numNodeLists);
  CHECK(riemannDvDx.size() == numNodeLists);

  // Derivative FieldLists.
  const auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  const auto  DrhoDx = derivatives.fields(GSPHFieldNames::densityGradient, Vector::zero);
  auto  normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DvolDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::volume, 0.0);
  //auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  //auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DmDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::mass, 0.0);
  auto  DEDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + GSPHFieldNames::thermalEnergy, 0.0);
  auto  DpDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + GSPHFieldNames::momentum, Vector::zero);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto  XSPHHfield = derivatives.fields("XSPHHfield", SymTensor::zero);
  
  auto  newRiemannDpDx = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannPressureGradient,Vector::zero);
  auto  newRiemannDvDx = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannVelocityGradient,Tensor::zero);

  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto& pairDepsDt = derivatives.getAny(HydroFieldNames::pairWork, vector<Scalar>());
  auto& pairMassFlux = derivatives.getAny(GSPHFieldNames::pairMassFlux, vector<Scalar>());
  CHECK(DrhoDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DvolDt.size() == numNodeLists);
  //CHECK(DvDt.size() == numNodeLists);
  //CHECK(DepsDt.size() == numNodeLists);
  CHECK(DEDt.size() == numNodeLists);
  CHECK(DpDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(newRiemannDpDx.size() == numNodeLists);
  CHECK(newRiemannDvDx.size() == numNodeLists);

  if (compatibleEnergy){
    pairAccelerations.resize(npairs);
    pairDepsDt.resize(2*npairs);
    pairMassFlux.resize(npairs);
  }

  this->computeMCorrection(time,dt,dataBase,state,derivatives);

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar psii,psij, Wi, gWi, Wj, gWj, Pstar, rhostari, rhostarj;
    Vector gradPsii, gradPsij, Ai, Aj, vstar;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    //auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);
    auto XSPHHfield_thread = XSPHHfield.threadCopy(threadStack);
    //auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DvolDt_thread = DvolDt.threadCopy(threadStack);
    auto DmDt_thread = DmDt.threadCopy(threadStack);
    auto DEDt_thread = DEDt.threadCopy(threadStack);
    auto DpDt_thread = DpDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto newRiemannDpDx_thread = newRiemannDpDx.threadCopy(threadStack);
    auto newRiemannDvDx_thread = newRiemannDvDx.threadCopy(threadStack);
    auto XSPHDeltaV_thread =  XSPHDeltaV.threadCopy(threadStack);
    auto normalization_thread = normalization.threadCopy(threadStack);
    
#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      const auto& mi = mass(nodeListi, i);
      const auto& mj = mass(nodeListj, j);

      if( mi >tiny or mj > tiny){

      // Get the state for node i.
      const auto& ui = nodalVelocity(nodeListi,i);
      const auto& riemannDpDxi = riemannDpDx(nodeListi, i);
      const auto& riemannDvDxi = riemannDvDx(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& voli = volume(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      //CHECK(mi > 0.0);
      CHECK(voli > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& normi = normalization_thread(nodeListi,i);
      //auto& DvDti = DvDt_thread(nodeListi, i);
      //auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvolDti = DvolDt_thread(nodeListi,i);
      auto& DmDti = DmDt_thread(nodeListi, i);
      auto& DEDti = DEDt_thread(nodeListi, i);
      auto& DpDti = DpDt_thread(nodeListi, i);
      auto& newRiemannDpDxi = newRiemannDpDx_thread(nodeListi, i);
      auto& newRiemannDvDxi = newRiemannDvDx_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& XSPHHfieldi = XSPHHfield_thread(nodeListi,i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi,i);
      const auto& gradRhoi = DrhoDx(nodeListi, i);
      const auto& Mi = M(nodeListi,i);


      // Get the state for node j
      const auto& uj = nodalVelocity(nodeListj,j);
      const auto& riemannDpDxj = riemannDpDx(nodeListj, j);
      const auto& riemannDvDxj = riemannDvDx(nodeListj, j);
      const auto& rj = position(nodeListj, j);
      //const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& volj = volume(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      //CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(volj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& normj = normalization_thread(nodeListj,j);
      //auto& DvDtj = DvDt_thread(nodeListj, j);
      //auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvolDtj = DvolDt_thread(nodeListj,j);
      auto& DmDtj = DmDt_thread(nodeListj, j);
      auto& DEDtj = DEDt_thread(nodeListj, j);
      auto& DpDtj = DpDt_thread(nodeListj, j);
      auto& newRiemannDpDxj = newRiemannDpDx_thread(nodeListj,j);
      auto& newRiemannDvDxj = newRiemannDvDx_thread(nodeListj,j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment(nodeListj, j);
      auto& XSPHHfieldj = XSPHHfield_thread(nodeListj,j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj,j);
      const auto& gradRhoj = DrhoDx(nodeListj, j);
      const auto& Mj = M(nodeListj,j);

      // Node displacement.
      const auto rij = ri - rj;
      const auto rMagij = rij.magnitude();
      const auto rMagij2 = rij.magnitude2();
      const auto rhatij =rij.unitVector();
      const auto vij = vi - vj;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);


      // Symmetrized kernel weight and gradient.
      //------------------------------------------------------
      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;

      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;

      psii = voli*Wi;
      psij = volj*Wj;
      gradPsii =  voli * Mi.Transpose()*gradWi;
      gradPsij =  volj * Mj.Transpose()*gradWj;

      const auto Astar = voli*gradPsii + volj*gradPsij;

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.
      const auto rij2 = rij.magnitude2();
      const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
      weightedNeighborSumi += std::abs(gWi);
      weightedNeighborSumj += std::abs(gWj);
      massSecondMomenti += gradWi.magnitude2()*thpt;
      massSecondMomentj += gradWj.magnitude2()*thpt;
      //const auto rij2 = rij.magnitude2();
      //const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
      //weightedNeighborSumi += std::abs(gWi);
      //weightedNeighborSumj += std::abs(gWj);
      // massSecondMomenti -= voli*rij.selfdyad()*gWi/rMagij;//.magnitude2()*thpt;
      // massSecondMomentj -= volj*rij.selfdyad()*gWj/rMagij;//.magnitude2()*thpt;

      // Determine an effective pressure including a term to fight the tensile instability.
      //const auto fij = epsTensile*pow(Wi/(Hdeti*WnPerh), nTensile);
      const auto fij = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
      const auto Ri = fij*(Pi < 0.0 ? -Pi : 0.0);
      const auto Rj = fij*(Pj < 0.0 ? -Pj : 0.0);
      const auto Peffi = Pi + Ri;
      const auto Peffj = Pj + Rj;

      // Reimann Solver and Fluxes
      //------------------------------------------------------
      // we'll clean this up when we have a gradient 
      // implementation we're in love with
      auto gradPi = riemannDpDxi;
      auto gradPj = riemannDpDxj;
      auto gradVi = riemannDvDxi;
      auto gradVj = riemannDvDxj;
      if (gradType==GradientType::SPHSameTimeGradient or
          gradType==GradientType::SPHUncorrectedGradient){
        gradPi = newRiemannDpDx(nodeListi,i);
        gradPj = newRiemannDpDx(nodeListj,j);
        gradVi = newRiemannDvDx(nodeListi,i);
        gradVj = newRiemannDvDx(nodeListj,j);
      }
      // need grad rho and grad eps
      riemannSolver.interfaceState(ri,           rj, 
                                   Hi,           Hj, 
                                   rhoi,         rhoj, 
                                   ci,           cj, 
                                   Peffi,        Peffj, 
                                   vi,           vj, 
                                   gradRhoi,     gradRhoj,
                                   gradPi,       gradPj, 
                                   gradVi,       gradVj, 
                                   Pstar,     //output
                                   vstar,     //output
                                   rhostari,  //output
                                   rhostarj); //output

      const auto fluxSwitch = 1;
      const auto vframe = (ui+uj)*0.5;
      const auto vflux = vstar-vframe;
      const auto fluxTowardsNodei = vflux.dot(rhatij) > 0;
      const auto rhostar = (fluxTowardsNodei ? rhostarj : rhostari); // we'll need to fix these later
      const auto epsstar = (fluxTowardsNodei ? epsj : epsi);         // we'll need to fix these later

      const auto massFlux = fluxSwitch * rhostar * vflux.dot(Astar);
      const auto momentumFlux = massFlux * vstar;
      const auto energyFlux = massFlux * epsstar;

      // mass
      //------------------------------------------------------
      DmDti -= massFlux;
      DmDtj += massFlux;

      // momentum
      //------------------------------------------------------
      const auto deltaDvDt = Pstar*Astar + momentumFlux; 
      DpDti -= deltaDvDt;
      DpDtj += deltaDvDt;

      // energy
      //------------------------------------------------------
      const auto deltaDepsDti = Pstar*Astar.dot(vi-vstar) - energyFlux;
      const auto deltaDepsDtj = Pstar*Astar.dot(vstar-vj) + energyFlux;

      //DepsDti += deltaDepsDti/mi;
      //DepsDtj += deltaDepsDtj/mj;

      DEDti += deltaDepsDti;
      DEDtj += deltaDepsDtj;
     
      if(compatibleEnergy){
        pairMassFlux[kk] = massFlux;
        pairAccelerations[kk] = deltaDvDt;
        pairDepsDt[2*kk]   = deltaDepsDti;
        pairDepsDt[2*kk+1] = deltaDepsDtj;
      }

      // volume change based on nodal velocity
      //-----------------------------------------------------
      DvolDti -= (ui-uj).dot(gradPsii);
      DvolDtj -= (ui-uj).dot(gradPsij);

      // gradients
      //------------------------------------------------------
      const auto deltaDvDxi = 2.0*(vi-vstar).dyad(gradPsii);
      const auto deltaDvDxj = 2.0*(vstar-vj).dyad(gradPsij);

      // based on riemann soln
      DvDxi -= deltaDvDxi;
      DvDxj -= deltaDvDxj;

      // while we figure out what we want ...
      switch(gradType){ 
        case GradientType::RiemannGradient: // default grad based on riemann soln
          newRiemannDvDxi -= deltaDvDxi;
          newRiemannDvDxj -= deltaDvDxj;
          newRiemannDpDxi -= 2.0*(Pi-Pstar)*gradPsii;
          newRiemannDpDxj -= 2.0*(Pstar-Pj)*gradPsij;
          break;
        case GradientType::HydroAccelerationGradient: // based on hydro accel for DpDx
          newRiemannDvDxi -= deltaDvDxi;
          newRiemannDvDxj -= deltaDvDxj;
          newRiemannDpDxi += deltaDvDt/voli;
          newRiemannDpDxj -= deltaDvDt/volj;
          break;
        case GradientType::SPHGradient: // raw gradients
          newRiemannDvDxi -= (vi-vj).dyad(gradPsii);
          newRiemannDvDxj -= (vi-vj).dyad(gradPsij);
          newRiemannDpDxi -= (Pi-Pj)*gradPsii;
          newRiemannDpDxj -= (Pi-Pj)*gradPsij;
          break;
        case GradientType::MixedMethodGradient: // raw gradient for P riemann gradient for v
          newRiemannDvDxi -= deltaDvDxi;
          newRiemannDvDxj -= deltaDvDxj;
          newRiemannDpDxi -= (Pi-Pj)*gradPsii;
          newRiemannDpDxj -= (Pi-Pj)*gradPsij;
          break;       
        default:
          break;
          // do nada    
        }

      // XSPH
      //-----------------------------------------------------------
      if (xsph) {
        XSPHDeltaVi -= psii*(vi-vj);
        XSPHDeltaVj -= psij*(vj-vi);
      }

      //if (sameMatij and diffuseEnergy){
          //linearReconstruction(ri,rj,epsi,epsj,DepsDxi,DepsDxj,epsLineari,epsLinearj);
          //const auto cijEff = max(min(cij + (vi-vj).dot(rhatij), cij),0.0);
          //const auto diffusion =  epsDiffusionCoeff*cijEff*(Hi-Hj)*etaij.dot(gradPsii)/(etaMagij*etaMagij+tiny);
      //XSPHHfieldi += psii*(massSecondMomentj);//ci*(Hi-Hj)*etai.dot(gradPsii)/(etaMagi*etaMagi+tiny);
      //XSPHHfieldj += psij*(massSecondMomenti);//cj*(Hj-Hi)*etaj.dot(gradPsij)/(etaMagj*etaMagj+tiny);i
      normi += psii;
      normj += psij;
      } //if statement

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
    const auto  kernelExtent = nodeList.neighbor().kernelExtent();
    const auto ni = nodeList.numInternalNodes();

#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      //const auto& mi = mass(nodeListi, i);
      const auto& voli = volume(nodeListi,i);
      //const auto& vi = velocity(nodeListi, i);
      const auto& ui = nodalVelocity(nodeListi,i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      //CHECK(mi > 0.0);
      CHECK(voli > 0.0);
      CHECK(Hdeti > 0.0);

      auto& normi = normalization(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DvolDti = DvolDt(nodeListi, i);
      //auto& DvDti = DvDt(nodeListi, i);     // FIX THIS
      //auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& XSPHHfieldi = XSPHHfield(nodeListi, i);
      
      normi += voli*Hdeti*W0;

      XSPHHfieldi += voli*Hdeti*W0*massSecondMomenti;
      //XSPHHfieldi /= max(normi,tiny);
      XSPHHfieldi /= Dimension::rootnu(std::max(XSPHHfieldi.Determinant(),tiny));

      DvolDti *= voli;

      // If needed finish the total energy derivative.
      //if (totalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      //weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));

      // Determine the position evolution, based on whether we're doing XSPH or not.
      DxDti = ui;
      if (xsph){
         DxDti += XSPHDeltaVi/max(tiny, normi);
      } 

      if(true){
      // The H tensor evolution.
      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;
      DHDti = smoothingScale.smoothingScaleDerivative(Hi,
                                                      ri,
                                                      DvDxi,
                                                      hmin,
                                                      hmax,
                                                      hminratio,
                                                      nPerh);
      
      //DHDti = 0.3*std::min(std::max((Ngb - 16.0),-1.0),1.0)/16.0 * Hi / dt;
      Hideali = smoothingScale.newSmoothingScale(Hi,
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
      }else{
      const auto Ngb_target = (Dimension::nDim == 3 ? 32 :
                                (Dimension::nDim == 2 ? 16 :
                                                        4));
      const auto C = (Dimension::nDim == 3 ? 1.33333*3.1415 :
                      (Dimension::nDim == 2 ? 3.1415         :
                                              1.0));
      // const auto detMSM = massSecondMomenti.Determinant();
      // weightedNeighborSumi = detMSM;//abs(XSPHHfieldi.Determinant()/max(Hdeti,tiny));
      // //const auto weightH = std::min(100*std::abs(weightedNeighborSumi),1.0);
      // //const auto Heffi = Hi ;
      // //if(abs(detMSM) > 1e-10){
      //   massSecondMomenti /= Dimension::rootnu(detMSM);
      
        
        // control that shape poopsie
        const auto circlerFactor = 0.00;//std::max(std::min(1.0,500.0*weightedNeighborSumi),0.05);//min(weightH,0.3);
        const auto stretchFactor = 0.4;

        const auto Ngb = C /(Hdeti*voli) * pow(kernelExtent,Dimension::nDim);
        //const auto Hstretch  =  circlerFactor * Dimension::rootnu(Hdeti)*SymTensor::one +
        //const auto Hstretch  =  ((1.00-stretchFactor)* SymTensor::one +
        //                               stretchFactor * XSPHHfieldi)*Hi;
        const auto scaleFactor = (1.0+0.5*(Ngb - Ngb_target)/Ngb_target);
        Hideali = std::min(std::max(scaleFactor,0.8),1.2) * Hi;
        // scale to enforce hmin/hmax
        DHDti = 0.25*(Hideali-Hi)/dt;
       }


      const auto C2 = (Dimension::nDim == 3 ? 1.33333*3.1415 :
                      (Dimension::nDim == 2 ? 3.1415         :
                                              1.0));
       weightedNeighborSumi = C2 /(Hdeti*voli) * pow(kernelExtent,Dimension::nDim);
      // } else{
      //   const auto stretchFactor = 0.00;
      //   const auto circlerFactor = 0.3;
      //   const auto Ngb = C /(Hdeti*voli) * pow(kernelExtent,Dimension::nDim);
      //   const auto  Hstretch  =  circlerFactor * Dimension::rootnu(Hdeti)*SymTensor::one +
      //                         ((1.00-stretchFactor-circlerFactor)*SymTensor::one)*Hi;
      //   const auto scaleFactor = (1.0+0.5*(Ngb - Ngb_target)/Ngb_target);
      //   Hideali = std::min(std::max(scaleFactor,0.9),1.1) * Hstretch;
      //   DHDti = 0.25*(Hideali-Hi)/dt;
      // }
    } // nodes loop
  }   // nodeLists loop
}     // eval derivs method 

//Ngb = 3.1415/(Hdeti*sum(Wi))
//dNdh = -3.1415/voliHdeti^2  dHdetidh - 3.1415/HdetiVoli^2 dVolidh
// Ngb = C h ^nu sum (W)
// dNdh = -N/h
//Hdeti = h3
//3h^2
//------------------------------------------------------------------------------
// EvalDerivs subroutine for spatial derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVHydroBase<Dimension>::
computeMCorrection(const typename Dimension::Scalar /*time*/,
                   const typename Dimension::Scalar /*dt*/,
                   const DataBase<Dimension>& dataBase,
                   const State<Dimension>& state,
                         StateDerivatives<Dimension>& derivatives) const {

  const auto calcSpatialGradients =  (this->gradientType() == GradientType::SPHSameTimeGradient 
                                  or  this->gradientType() == GradientType::SPHUncorrectedGradient);
  const auto correctSpatialGradients = (this->gradientType() == GradientType::SPHSameTimeGradient);
  // The kernels and such.
  const auto& W = this->kernel();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists. 
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);

  CHECK(massDensity.size() == numNodeLists);
  CHECK(volume.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);

  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  DrhoDx = derivatives.fields(GSPHFieldNames::densityGradient, Vector::zero);
  auto  newRiemannDpDx = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannPressureGradient,Vector::zero);
  auto  newRiemannDvDx = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + GSPHFieldNames::RiemannVelocityGradient,Tensor::zero);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  
  CHECK(M.size() == numNodeLists);
  CHECK(DrhoDx.size() == numNodeLists);
  CHECK(newRiemannDpDx.size() == numNodeLists);
  CHECK(newRiemannDvDx.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto M_thread = M.threadCopy(threadStack);
    auto DrhoDx_thread = DrhoDx.threadCopy(threadStack);
    auto newRiemannDpDx_thread = newRiemannDpDx.threadCopy(threadStack);
    auto newRiemannDvDx_thread = newRiemannDvDx.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);
#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;
      
      // Get the state for node i.
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& ri = position(nodeListi, i);
      const auto& voli = volume(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(voli > 0.0);
      CHECK(Hdeti > 0.0);

      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);
      auto& DrhoDxi = DrhoDx_thread(nodeListi, i);
      auto& Mi = M_thread(nodeListi, i);

      // Get the state for node j
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& rj = position(nodeListj, j);
      const auto& volj = volume(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(volj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);
      auto& DrhoDxj = DrhoDx_thread(nodeListj, j);
      auto& Mj = M_thread(nodeListj, j);

      const auto rij = ri - rj;
      const auto rMagij = safeInv(rij.magnitude());

      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      const auto gWi = W.gradValue(etaMagi, Hdeti);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;

      const auto gWj = W.gradValue(etaMagj, Hdetj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;

      const auto gradPsii = voli*gradWi;
      const auto gradPsij = volj*gradWj;

      Mi -= rij.dyad(gradPsii);
      Mj -= rij.dyad(gradPsij);
      
      DrhoDxi -= (rhoi - rhoj) * gradPsii;
      DrhoDxj -= (rhoi - rhoj) * gradPsij;
      
      //massSecondMomenti -= voli*rij.selfdyad()*gWi*rMagij;//.magnitude2()*thpt;
      //massSecondMomentj -= volj*rij.selfdyad()*gWj*rMagij;//.magnitude2()*thpt;

      // // based on nodal values
      if (calcSpatialGradients){
        const auto& vi = velocity(nodeListi, i);
        const auto& Pi = pressure(nodeListi, i);
        const auto& vj = velocity(nodeListj, j);
        const auto& Pj = pressure(nodeListj, j);
        auto& newRiemannDpDxi = newRiemannDpDx_thread(nodeListi, i);
        auto& newRiemannDvDxi = newRiemannDvDx_thread(nodeListi, i);
        auto& newRiemannDpDxj = newRiemannDpDx_thread(nodeListj, j);
        auto& newRiemannDvDxj = newRiemannDvDx_thread(nodeListj, j);

        newRiemannDpDxi -= (Pi-Pj)*gradPsii;
        newRiemannDpDxj -= (Pi-Pj)*gradPsij;

        newRiemannDvDxi -= (vi-vj).dyad(gradPsii);
        newRiemannDvDxj -= (vi-vj).dyad(gradPsij);
      }
    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  
  // Finish up the spatial gradient calculation
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = M[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      //auto& massSecondMomenti = massSecondMoment(nodeListi,i);
      
      //const auto detMSM = massSecondMomenti.Determinant();
      //massSecondMomenti /= Dimension::rootnu(detMSM);

      const auto Mdeti = std::abs(Mi.Determinant());

      const auto enoughNeighbors =  numNeighborsi > Dimension::pownu(2);
      const auto goodM =  (Mdeti > 1e-2 and enoughNeighbors);                   

      Mi = ( goodM ? Mi.Inverse() : Tensor::one);

      if (correctSpatialGradients){
        auto& newRiemannDpDxi = newRiemannDpDx(nodeListi, i);
        auto& newRiemannDvDxi = newRiemannDvDx(nodeListi, i);

        newRiemannDpDxi = Mi.Transpose()*newRiemannDpDxi;
        newRiemannDvDxi = newRiemannDvDxi*Mi;
      }
    }
    
  }
  
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
         boundItr != this->boundaryEnd();
         ++boundItr)(*boundItr)->applyFieldListGhostBoundary(M);

  if (calcSpatialGradients){ 
    for (ConstBoundaryIterator boundItr = this->boundaryBegin();
          boundItr != this->boundaryEnd();
           ++boundItr){
      (*boundItr)->applyFieldListGhostBoundary(newRiemannDpDx);
      (*boundItr)->applyFieldListGhostBoundary(newRiemannDvDx);
    }
  }
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
 
}

} // spheral namespace
