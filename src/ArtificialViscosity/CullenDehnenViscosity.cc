//---------------------------------Spheral++----------------------------------//
// Cullen adn Dehnen Viscosity 
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DataOutput/Restart.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "IncrementCullenMultipliers.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Material/EquationOfState.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"

#include "CullenDehnenViscosity.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {
    
//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
CullenDehnenViscosity<Dimension>::
CullenDehnenViscosity(ArtificialViscosity<Dimension>& q,
                      const TableKernel<Dimension>& W,
                      const Scalar alphMax, //Parameter = 2.0 in Hopkins 2014 and Price 2004, = 1.5 in Rosswog 2000
                      const Scalar alphMin, //Parameter = 0.02 in Hopkins 2014  
                      const Scalar betaC, //Parameter = 0.7 Hopkins 2014
                      const Scalar betaD, //Parameter = 0.05 Hopkins 2014
                      const Scalar betaE, //Parameter = 1.0 in Hopkins 2014, = 2.0 in Cullen 2010
                      const Scalar fKern, //Parameter = 1/3 Hopkins 2014 for quinitc spline
                      const bool boolHopkins): //use Hopkins Reformulation
  Physics<Dimension>(),
  mPrevDvDt(FieldStorageType::CopyFields),
  mPrevDivV(FieldStorageType::CopyFields),
  mCullAlpha(FieldStorageType::CopyFields),
  mPrevDivV2(FieldStorageType::CopyFields),
  mCullAlpha2(FieldStorageType::CopyFields),
  mDalphaDt(FieldStorageType::CopyFields),
  mAlphaLocal(FieldStorageType::CopyFields),
  malphMax(alphMax),
  malphMin(alphMin),
  mbetaC(betaC),
  mbetaD(betaD),
  mbetaE(betaE),
  mfKern(fKern),
  mboolHopkins(boolHopkins),
  myq(q),
  mKernel(W),
  mRestart(registerWithRestart(*this)){
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CullenDehnenViscosity<Dimension>::
~CullenDehnenViscosity() {
}

    
// Accessor Fns
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
alphMax(Scalar val)
{
    malphMax = val;
}

template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
alphMin(Scalar val)
{
    malphMin = val;
}

template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
betaE(Scalar val)
{
    mbetaE = val;
}

template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
betaD(Scalar val)
{
    mbetaD = val;
}

template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
betaC(Scalar val)
{
    mbetaC = val;
}

template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
fKern(Scalar val)
{
    mfKern = val;
}

template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
boolHopkins(bool val)
{
    mboolHopkins = val;
}

//------------------------------------------------------------------------------
// Access the main kernel
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
CullenDehnenViscosity<Dimension>::kernel() const {
  return mKernel;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CullenDehnenViscosity<Dimension>::PrevDvDt() const {
   return mPrevDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::PrevDivV() const {
   return mPrevDivV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::PrevDivV2() const {
   return mPrevDivV2;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::CullAlpha() const {
   return mCullAlpha;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::CullAlpha2() const {
   return mCullAlpha2;
}

template<typename Dimension>
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
alphMax() const{return malphMax;}

template<typename Dimension>
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
alphMin() const{return malphMin;}

template<typename Dimension>
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
betaE() const{return mbetaE;}

template<typename Dimension>
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
betaD() const{return mbetaD;}

template<typename Dimension>
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
betaC() const{return mbetaC;}

template<typename Dimension>
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
fKern() const{return mfKern;}

template<typename Dimension>
bool CullenDehnenViscosity<Dimension>::
boolHopkins() const{return mboolHopkins;}

template<typename Dimension>
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::
DalphaDt() const{ return mDalphaDt;}

template<typename Dimension>
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::
alphaLocal() const{ return mAlphaLocal;}
    
//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  mPrevDvDt = dataBase.newFluidFieldList(Vector::zero, "mPrevDvDt");
  mPrevDivV = dataBase.newFluidFieldList(0.0, "mPrevDivV");
  mCullAlpha = dataBase.newFluidFieldList(1.0, "mCullAlpha");
  mPrevDivV2 = dataBase.newFluidFieldList(0.0, "mPrevDivV2");
  mCullAlpha2 = dataBase.newFluidFieldList(1.0, "mCullAlpha2");
  mDalphaDt = dataBase.newFluidFieldList(0.0, "Cullen alpha delta");
  mAlphaLocal = dataBase.newFluidFieldList(0.0, "Cullen alpha local");

  FieldList<Dimension, Scalar>& rvQ = myq.CqMultiplier();
  FieldList<Dimension, Scalar>& rvL = myq.ClMultiplier();
  rvQ = dataBase.newFluidFieldList(malphMin, HydroFieldNames::ArtificialViscousCqMultiplier);  // This will override the Q initializer intializing these to unity.
  rvL = dataBase.newFluidFieldList(malphMin, HydroFieldNames::ArtificialViscousClMultiplier);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;
  dataBase.resizeFluidFieldList(mPrevDvDt, Vector::zero, "mPrevDvDt", false);
  dataBase.resizeFluidFieldList(mPrevDivV, 0.0, "mPrevDivV", false);
  dataBase.resizeFluidFieldList(mCullAlpha, 1.0, "mCullAlpha", false);
  state.enroll(mPrevDvDt);
  state.enroll(mPrevDivV);
  state.enroll(mCullAlpha);

  FieldList<Dimension, Scalar>& rvAlphaQ = myq.CqMultiplier();
  FieldList<Dimension, Scalar>& rvAlphaL = myq.ClMultiplier();
  PolicyPointer alphaPolicy(new IncrementCullenMultipliers<Dimension>(malphMin, malphMax, mboolHopkins));
  state.enroll(rvAlphaQ, alphaPolicy);
  state.enroll(rvAlphaL);
}
    
//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  dataBase.resizeFluidFieldList(mPrevDivV2, 0.0, "mPrevDivV2", false);
  dataBase.resizeFluidFieldList(mCullAlpha2, 1.0, "mCullAlpha2", false);
  derivs.enroll(mPrevDivV2);
  derivs.enroll(mCullAlpha2);
  derivs.enroll(mDalphaDt);
  derivs.enroll(mAlphaLocal);
}

template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
}

//------------------------------------------------------------------------------
// Determine the Cullen Coefficients
//------------------------------------------------------------------------------
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
finalizeDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  // The kernels
  const TableKernel<Dimension>& W = this->kernel();
  const Scalar kernelExtent = W.kernelExtent();

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();
    
  //State Fluid Lists
  const FieldList<Dimension, Scalar> reducingViscosityMultiplierQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);

  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  
  FieldList<Dimension, Vector> prevDvDt = state.fields("mPrevDvDt", Vector::zero);

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
    
  CHECK(reducingViscosityMultiplierQ.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> alpha_local = derivs.fields("Cullen alpha local", 0.0);
  FieldList<Dimension, Scalar> DalphaDt = derivs.fields("Cullen alpha delta", 0.0);

  // We're using the hydro derivatives.
  const FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);

  // Apply boundaries to DvDx.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(DvDx);
  }
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  // We need to compute the R factor, which involves walking the neighbors.  We can simultaneously compute the signal velocity.
  FieldList<Dimension, Scalar> R = dataBase.newFluidFieldList(0.0, "Cullen R limiter");
  FieldList<Dimension, Scalar> vsig = dataBase.newFluidFieldList(0.0, "Cullen signal velocity");
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const int firstGhostNodei = DvDx[nodeListi]->nodeList().firstGhostNode();

    // Iterate over the nodes in this node list.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar mi = mass(nodeListi, i);
      const Scalar rhoi = massDensity(nodeListi, i);
      const Scalar divvi = DvDx(nodeListi, i).Trace();
      const Scalar csi = soundSpeed(nodeListi, i);
      Scalar& Ri = R(nodeListi, i);
      Scalar& vsigi = vsig(nodeListi, i);

      // Neighbors!
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Walk the neighbor nodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      
        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = DvDx[nodeListj]->nodeList().firstGhostNode();

          // Loop over the neighbors.
#if defined __INTEL_COMPILER
#pragma vector always
#endif
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            
            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {

              // Get the state for node i.
              const Vector& rj = position(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar mj = mass(nodeListj, j);
              const Scalar divvj = DvDx(nodeListj, j).Trace();
              const Scalar csj = soundSpeed(nodeListj, j);
              Scalar& Rj = R(nodeListj, j);
              Scalar& vsigj = vsig(nodeListj, j);

              // Symmetrized kernel weight and gradient.
              const Vector rij = ri - rj;
              const Scalar Wi = W(Hi*rij, Hdeti);
              const Scalar Wj = W(Hj*rij, Hdetj);

              // Compute R.
              Ri += sgn(divvj)*mj*Wi;
              Rj += sgn(divvi)*mi*Wj;

              // Check the signal velocity.
              const Scalar vsigij = 0.5*(csi + csj) - std::min(0.0, (vi - vj).dot(rij.unitVector()));
              vsigi = std::max(vsigi, vsigij);
              vsigj = std::max(vsigj, vsigij);
            }
          }
        }
      }

      // Finish Ri.
      Ri /= rhoi;
    }
  }

  // Now we have DvDx, DvDtDx, R, and vsig.  We can compute the Cullen & Dehnen viscosity coefficients.
  // We also repurpose 
  //   mCullAlpha  -> alpha0
  //   mCullAlpha2 -> alpha_tmp.
  //   alpha_local -> alpha
  const FieldList<Dimension, Scalar> prevDivV = state.fields("mPrevDivV", 0.0);
  const FieldList<Dimension, Scalar> alpha0 = state.fields("mCullAlpha", 0.0);
  FieldList<Dimension, Scalar> alpha_tmp = derivs.fields("mCullAlpha2", 0.0);
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const size_t n = DvDx[nodeListi]->numInternalElements();
    for (size_t i = 0; i != n; ++i) {
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Ri = R(nodeListi, i);
      const Scalar vsigi = vsig(nodeListi, i);
      const Tensor DvDxi = DvDx(nodeListi, i);
      const Scalar divvi = DvDxi.Trace();
      const Scalar divai = (divvi - prevDivV(nodeListi, i))*safeInvVar(dt);
      const Tensor Si = DvDxi.Symmetric() - divvi/Dimension::nDim * Tensor::one;
      const Scalar alphai = reducingViscosityMultiplierQ(nodeListi, i);
      Scalar& alpha_locali = alpha_local(nodeListi, i);
      if (mboolHopkins) {
        const Scalar hi = kernelExtent*Dimension::nDim/Hi.Trace();  // Harmonic averaging.
        const Scalar alpha0i = alpha0(nodeListi, i);
        Scalar& alpha_tmpi = alpha_tmp(nodeListi, i);
        alpha_tmpi = ((divvi >= 0.0 or divai >= 0.0) ?
                      0.0 :
                      malphMax*std::abs(divai)*safeInvVar(std::abs(divai) + mbetaC*vsigi*vsigi/FastMath::pow2(mfKern*hi)));
        const Scalar thpt = FastMath::pow2(mbetaE*FastMath::pow4(1.0 - Ri)*divvi);
        // const Scalar taui = hi/kernelExtent*safeInvVar(2.0*mbetaD*vsigi);
        alpha_locali = std::max(malphMin, thpt*alpha0i*safeInvVar(thpt + (Si*Si.Transpose()).Trace()));
        // We abuse DalphaDt here to store the new value for alpha0 in the Hopkins approximation.
        DalphaDt(nodeListi, i) = alpha_tmpi + std::max(0.0, alpha0i - alpha_tmpi)*exp(-mbetaD*dt*abs(vsigi)/(2.0*mfKern*hi));
      } else {
        const Scalar hi = kernelExtent*Dimension::nDim/Hi.Trace();  // Harmonic averaging.
        const Scalar thpt = FastMath::pow2(2.0*FastMath::pow4(1.0 - Ri)*divvi);
        const Scalar zetai = thpt*safeInvVar(thpt + (Si*Si.Transpose()).Trace());
        const Scalar Ai = zetai*std::max(-divai, 0.0);
        const Scalar taui = hi*safeInvVar(2.0*mbetaD*vsigi);
        DalphaDt(nodeListi, i) = std::min(0.0, alpha_locali - alphai)*safeInv(taui);
        alpha_locali = std::max(alphai, malphMax*hi*hi*Ai*safeInvVar(vsigi*vsigi + hi*hi*Ai));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Finalize
//------------------------------------------------------------------------------
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {
  FieldList<Dimension, Vector> prevDvDt = state.fields("mPrevDvDt", Vector::zero);
  const FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  prevDvDt.assignFields(DvDt);
  FieldList<Dimension, Scalar> prevDivV = state.fields("mPrevDivV", 0.0);
  const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  const unsigned numNodeLists = DvDx.numFields();
  CHECK(prevDivV.numFields() == numNodeLists);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = DvDx[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      prevDivV(nodeListi, i) = DvDx(nodeListi, i).Trace();
    }
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
    file.write(mPrevDvDt, pathName + "/prevDvDt");
    file.write(mPrevDivV, pathName + "/prevDivV");
    file.write(mCullAlpha, pathName + "/cullAlpha");
    file.write(mPrevDivV2, pathName + "/prevDivV2");
    file.write(mCullAlpha2, pathName + "/cullAlpha2");
    file.write(mDalphaDt, pathName + "/DalphaDt");
    file.write(mAlphaLocal, pathName + "/alphaLocal");
}
    
//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
    file.read(mPrevDvDt, pathName + "/prevDvDt");
    file.read(mPrevDivV, pathName + "/prevDivV");
    file.read(mCullAlpha, pathName + "/cullAlpha");
    file.read(mPrevDivV2, pathName + "/prevDivV2");
    file.read(mCullAlpha2, pathName + "/cullAlpha2");
    file.read(mDalphaDt, pathName + "/DalphaDt");
    file.read(mAlphaLocal, pathName + "/alphaLocal");
}
    
//------------------------------------------------------------------------------
// Boundaries
//------------------------------------------------------------------------------
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs)
{
    FieldList<Dimension, Vector> prevDvDt = state.fields("mPrevDvDt", Vector::zero);
    FieldList<Dimension, Scalar> prevDivV = state.fields("mPrevDivV", 0.0);
    FieldList<Dimension, Scalar> cullAlpha = state.fields("mCullAlpha", 0.0);
    FieldList<Dimension, Scalar> prevDivV2 = derivs.fields("mPrevDivV2", 0.0);
    FieldList<Dimension, Scalar> cullAlpha2 = derivs.fields("mCullAlpha2", 0.0);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
      (*boundaryItr)->applyFieldListGhostBoundary(prevDvDt);
      (*boundaryItr)->applyFieldListGhostBoundary(prevDivV);
      (*boundaryItr)->applyFieldListGhostBoundary(cullAlpha);
      (*boundaryItr)->applyFieldListGhostBoundary(prevDivV2);
      (*boundaryItr)->applyFieldListGhostBoundary(cullAlpha2);
    }

}
    
template<typename Dimension>
void
CullenDehnenViscosity<Dimension>::
enforceBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs)
{
    FieldList<Dimension, Vector> prevDvDt = state.fields("mPrevDvDt", Vector::zero);
    FieldList<Dimension, Scalar> prevDivV = state.fields("mPrevDivV", 0.0);
    FieldList<Dimension, Scalar> cullAlpha = state.fields("mCullAlpha", 0.0);
    FieldList<Dimension, Scalar> prevDivV2 = derivs.fields("mPrevDivV2", 0.0);
    FieldList<Dimension, Scalar> cullAlpha2 = derivs.fields("mCullAlpha2", 0.0);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
      (*boundaryItr)->enforceFieldListBoundary(prevDvDt);
      (*boundaryItr)->enforceFieldListBoundary(prevDivV);
      (*boundaryItr)->enforceFieldListBoundary(cullAlpha);
      (*boundaryItr)->enforceFieldListBoundary(prevDivV2);
      (*boundaryItr)->enforceFieldListBoundary(cullAlpha2);
    }
}
    
//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename CullenDehnenViscosity<Dimension>::TimeStepType
CullenDehnenViscosity<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
    return TimeStepType(1.0e100, "Rate of viscosity change -- NO VOTE.");
}

}
