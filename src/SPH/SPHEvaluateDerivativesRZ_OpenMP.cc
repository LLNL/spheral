namespace Spheral {

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
  const double rhoTiny = 1.0e-10;
  const double W0 = W(0.0, 1.0);

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
  FieldList<Dimension, Vector> DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
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
    const auto ni = connectivityMap.numNodes(nodeListi);
#pragma omp parallel for
    for (auto iItr = 0; iItr < ni; ++iItr) {
      const auto i = connectivityMap.ithNode(nodeListi, iItr);

      // Prepare to accumulate the time.
      const Time start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const Vector& posi = position(nodeListi, i);
      const Scalar ri = abs(posi.y());
      const Scalar circi = 2.0*M_PI*ri;
      const Scalar mi = mass(nodeListi, i);
      const Scalar mRZi = mi/circi;
      const Vector& vi = velocity(nodeListi, i);
      const Scalar rhoi = massDensity(nodeListi, i);
      const Scalar epsi = specificThermalEnergy(nodeListi, i);
      const Scalar Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar ci = soundSpeed(nodeListi, i);
      const Scalar& omegai = omega(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar safeOmegai = safeInv(omegai, tiny);
      const Scalar zetai = abs((Hi*posi).y());
      const Scalar hri = ri*safeInv(zetai);
      const Scalar riInv = safeInv(ri, 0.25*hri);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

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
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Get the state for node j
            const Vector& posj = position(nodeListj, j);
            const Scalar rj = abs(posj.y());
            const Scalar circj = 2.0*M_PI*rj;
            const Scalar mj = mass(nodeListj, j);
            const Scalar mRZj = mj/circj;
            const Vector& vj = velocity(nodeListj, j);
            const Scalar rhoj = massDensity(nodeListj, j);
            const Scalar epsj = specificThermalEnergy(nodeListj, j);
            const Scalar Pj = pressure(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar cj = soundSpeed(nodeListj, j);
            const Scalar& omegaj = omega(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();
            const Scalar safeOmegaj = safeInv(omegaj, tiny);
            const Scalar zetaj = abs((Hj*posj).y());
            CHECK(rhoj > 0.0);
            CHECK(Hdetj > 0.0);

            // Flag if this is a contiguous material pair or not.
            const bool sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

            // Node displacement.
            const Vector xij = posi - posj;
            const Vector etai = Hi*xij;
            const Vector etaj = Hj*xij;
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
            const double fweightij = sameMatij ? 1.0 : mRZj*rhoi/(mRZi*rhoj);
            const double xij2 = xij.magnitude2();
            const auto thpt = xij.selfdyad()*safeInvVar(xij2*xij2*xij2);
            weightedNeighborSumi +=     fweightij*std::abs(gWi);
            massSecondMomenti +=     fweightij*gradWi.magnitude2()*thpt;

            // Contribution to the sum density.
            if (nodeListi == nodeListj) {
              rhoSumi += mRZj*Wi;
              normi += mRZi/rhoi*Wi;
            }

            // Compute the pair-wise artificial viscosity.
            const Vector vij = vi - vj;
            const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                      posi, etai, vi, rhoi, ci, Hi,
                                                      posj, etaj, vj, rhoj, cj, Hj);
            const Vector Qacci = 0.5*(QPiij.first *gradWQi);
            const Vector Qaccj = 0.5*(QPiij.second*gradWQj);
            const Scalar workQi = 0.5*(QPiij.first *vij).dot(gradWQi);
            // const Scalar workQi = vij.dot(Qacci);
            const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
            maxViscousPressurei = max(maxViscousPressurei, Qi);
            effViscousPressurei += mRZj/rhoj * Qi * Wi;
            viscousWorki += mRZj*workQi;

            // Acceleration.
            CHECK(rhoi > 0.0);
            CHECK(rhoj > 0.0);
            const double Prhoi = safeOmegai*Pi/(rhoi*rhoi);
            const double Prhoj = safeOmegaj*Pj/(rhoj*rhoj);
            const Vector deltaDvDt = (Prhoi*gradWi + Prhoj*gradWj) + Qacci + Qaccj;
            DvDti -= mRZj*deltaDvDt;

            // Specific thermal energy evolution.
            DepsDti += mRZj*(Prhoi*vij.dot(gradWi) + workQi);
            if (mCompatibleEnergyEvolution) {
              pairAccelerationsi.push_back(-mRZj*deltaDvDt);
              pairAccelerationsi.push_back( mRZi*deltaDvDt);
            }

            // Velocity gradient.
            const Tensor deltaDvDxi = mRZj*vij.dyad(gradWi);
            DvDxi -= deltaDvDxi;
            if (sameMatij) {
              localDvDxi -= deltaDvDxi;
            }

            // Estimate of delta v (for XSPH).
            if (sameMatij or min(zetai, zetaj) < 1.0) {
              const double wXSPHij = 0.5*(mRZi/rhoi*Wi + mRZj/rhoj*Wj);
              XSPHWeightSumi += wXSPHij;
              XSPHDeltaVi -= wXSPHij*vij;
            }

            // Linear gradient correction term.
            Mi -= mRZj*xij.dyad(gradWi);
            if (sameMatij) {
              localMi -= mRZj*xij.dyad(gradWi);
            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not mCompatibleEnergyEvolution or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == 2*numNeighborsi));

      // Get the time for pairwise interactions.
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/(numNeighborsi + 1.0e-30);

      // Add the self-contribution to density sum.
      rhoSumi += mRZi*W0*Hdeti;
      rhoSumi /= circi;
      normi += mRZi/rhoi*W0*Hdeti;

      // Finish the acceleration.
      pairAccelerationsi.push_back(Vector::zero);

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

      // Finish the continuity equation.
      XSPHWeightSumi += Hdeti*mRZi/rhoi*W0;
      CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
      XSPHDeltaVi /= XSPHWeightSumi;
      const Scalar vri = vi.y(); // + XSPHDeltaVi.y();
      DrhoDti = -rhoi*(DvDxi.Trace() + vri*riInv);

      // Finish the specific thermal energy evolution.
      DepsDti -= Pi/rhoi*vri*riInv;

      // If needed finish the total energy derivative.
      if (mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             posi,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
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

      // Increment the work for i.
      worki += Timing::difference(start, Timing::currentTime());
    }
  }
}

}
