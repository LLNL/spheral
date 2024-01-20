/*
  HelmholtzEquationOfState.cc

  Created by Raskin, Cody Dantes on 8/28/14.

  The helmholtz EOS solves the helmholtz free energy equation by means of
  a table lookup. This table file "helm_table.dat" MUST be present at runtime
  in the same location as the run script or the EOS will fail.
 
*/

#include "HelmholtzEquationOfState.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include <memory>

#define Fortran2(x) x##_
#define NBLOCK 100

extern "C" {
  void Fortran2(init_helm_table)();
        
  void Fortran2(get_helm_table)(double *f,double *fd,double *ft,double *fdd,
                                double *ftt,double *fdt,double *fddt,
                                double *fdtt,double *fddtt,double *dpdf,
                                double *dpdfd,double *dpdft,double *dpdfdt,
                                double *ef,double *efd,double *eft,double *efdt,
                                double *xf,double *xfd,double *xft,double *xfdt);
        
  void Fortran2(wrapper_invert_helm_ed)(int *npart, double *density,
                                        double *energy, double *abar,
                                        double *zbar, double *temperature,
                                        double *pressure, double *small_temp, double *vsound,
                                        double *gamma, double *entropy);
        
  void Fortran2(wrapper_helmeos)(int *npart, double *den_row,
                                 double *etot_row, double *abar_row,
                                 double *zbar_row, double *temperature,
                                 double *pressure);
        
  void Fortran2(set_helm_table)(double *f, double *fd, double *ft, double *fdd,
                                double *ftt, double *fdt, double *fddt,
                                double *fdtt, double *fddtt, double *dpdf,
                                double *dpdfd, double *dpdft, double *dpdfdt,
                                double *ef, double *efd, double *eft,
                                double *efdt, double *xf, double *xfd,
                                double *xft, double *xfdt);
        
  void Fortran2(azbar)(double *xmass, double *aion, double *zion, int *ionmax,
                       double *ymass, double *abar, double *zbar);
}

namespace Spheral {

using std::shared_ptr;

template<typename Dimension>
HelmholtzEquationOfState<Dimension>::
HelmholtzEquationOfState(const PhysicalConstants& constants,
                         const double minimumPressure,
                         const double maximumPressure,
                         const double minimumTemperature,
                         const MaterialPressureMinType minPressureType,
                         const Scalar abar0,
                         const Scalar zbar0,
                         const double externalPressure):
  EquationOfState<Dimension>(constants, minimumPressure, maximumPressure, minPressureType, externalPressure),
  myAbar(),
  myZbar(),
  mySpecificThermalEnergy(),
  myMassDensity(),
  myTemperature(),
  myPressure(),
  mySoundSpeed(),
  myGamma(),
  myEntropy(),
  mabar0(abar0),
  mzbar0(zbar0),
  mPmin(minimumPressure),
  mPmax(maximumPressure),
  mTmin(minimumTemperature),
  mConstants(constants)
{
  needUpdate = 1; // flip this on and off later
  //Fortran2(init_helm_table);
      
  /*
    need to set constants here to use CGS no matter what
    is passed in
  */
      
  mDistincm           = mConstants.unitLengthMeters() / 0.01;
  mMassing            = mConstants.unitMassKg() / 0.001;
  mDensingpccm        = mMassing / (mDistincm*mDistincm*mDistincm);
  mTimeins            = mConstants.unitTimeSec();
  mVelincmps          = mDistincm / mTimeins;
  mEnergyinergpg      = mVelincmps*mVelincmps;
  mPressureinbarye    = mMassing / (mDistincm*mTimeins*mTimeins);
      
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
HelmholtzEquationOfState<Dimension>::
~HelmholtzEquationOfState() {
}


/* check mzbar != nullptr?? */


//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
      
  storeFields(massDensity, specificThermalEnergy);
      
  /* what follows is a method of breaking up the input to the fortran routine into smaller chunks */
      
  int npart = massDensity.numElements();
  int nblock = NBLOCK;
  int nrest = npart % nblock;
  int nloop = (npart - nrest)/nblock;
  unsigned int k = 0;
      
  if(needUpdate){
    for (int j=0; j<nloop; ++j)
    {
      k = j*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nblock, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
    /* now do the rest */
    if (nrest > 0)
    {
      k = nloop*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nrest, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
          
  }
      
  for (int i = 0; i != npart; ++i) {
    Pressure(i) = (*myPressure)[i] / mPressureinbarye;
  }
}

//------------------------------------------------------------------------------
// Set the pressure and derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& Pressure,
                     Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                     Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {
  VERIFY2(false, "HelmholtzEquationOfState::setPressureAndDerivatives not supported");
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());

  storeFields(massDensity, specificThermalEnergy);
      
  /* what follows is a method of breaking up the input to the fortran routine into smaller chunks */
      
  int npart = massDensity.numElements();
  int nblock = NBLOCK;
  int nrest = npart % nblock;
  int nloop = (npart - nrest)/nblock;
  unsigned int k = 0;
      
  if(needUpdate){
    for (int j=0; j<nloop; ++j)
    {
      k = j*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nblock, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
    /* now do the rest */
    if (nrest > 0)
    {
      k = nloop*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nrest, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
          
  }

  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    temperature(i) = (*myTemperature)[i];
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
  
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& /*specificThermalEnergy*/,
                         const Field<Dimension, Scalar>& /*massDensity*/,
                         const Field<Dimension, Scalar>& /*temperature*/) const {
  /*
      
    CHECK(valid());
      
    storeFields(massDensity, specificThermalEnergy);
      
    int npart       = massDensity.numElements();
    myMassDensity   = shared_ptr<Field<Dimension, Scalar> >(&massDensity);
    myTemperature   = shared_ptr<Field<Dimension, Scalar> >(&temperature);
      
    if(needUpdate){
    Fortran2(wrapper_helmeos)(&npart, &(myMassDensity->at(0)), &(mySpecificThermalEnergy->at(0)),
    &(myAbar->at(0)), &(myZbar->at(0)), &(myTemperature->at(0)),
    &(myPressure->at(0)));
    }
      
    for (size_t i = 0; i != npart; ++i) {
    specificThermalEnergy(i) = mySpecificThermalEnergy->at(i);
    }
  */
      
  VERIFY2(false,"Helmhotlz EOS does not support this call.");
}
  

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& /*specificHeat*/,
                const Field<Dimension, Scalar>& /*massDensity*/,
                const Field<Dimension, Scalar>& /*temperature*/) const {
  /*  
      CHECK(valid());
      
      storeFields(massDensity, specificThermalEnergy);
      
      const double kB = mConstants.kB();
      const double mp = mConstants.protonMass();
      int npart = myGamma->numElements();
      double Cv;
      
      for (size_t i = 0; i != npart; ++i)
      Cv += kB/(myGamma->at(i)*myAbar->at(i)*mp);
      specificHeat = Cv/npart;

  */
  VERIFY2(false,"Helmholtz EOS does not support this call.");
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());

  storeFields(massDensity, specificThermalEnergy);
      
  /* what follows is a method of breaking up the input to the fortran routine into smaller chunks */
      
  int npart = massDensity.numElements();
  int nblock = NBLOCK;
  int nrest = npart % nblock;
  int nloop = (npart - nrest)/nblock;
  unsigned int k = 0;
      
  if(needUpdate){
    for (int j=0; j<nloop; ++j)
    {
      k = j*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nblock, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
    /* now do the rest */
    if (nrest > 0)
    {
      k = nloop*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nrest, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
          
  }
      
  for (int i = 0; i != npart; ++i) {
    soundSpeed(i) = (*mySoundSpeed)[i] / mVelincmps;
    (*myGamma)[i] = soundSpeed(i) * soundSpeed(i) * massDensity(i) / ((*myPressure)[i] / mPressureinbarye);
  }
}
  
//------------------------------------------------------------------------------
// Set gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
      
  storeFields(massDensity, specificThermalEnergy);
      
  /* what follows is a method of breaking up the input to the fortran routine into smaller chunks */
      
  int npart = massDensity.numElements();
  int nblock = NBLOCK;
  int nrest = npart % nblock;
  int nloop = (npart - nrest)/nblock;
  unsigned int k = 0;
      
  if(needUpdate){
    for (int j=0; j<nloop; ++j)
    {
      k = j*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nblock, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
    /* now do the rest */
    if (nrest > 0)
    {
      k = nloop*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nrest, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
          
  }
      
  for (int i = 0; i != npart; ++i) {
    gamma(i) = (*myGamma)[i];
  }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for a
// Helmholtz gas.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
      
  setPressure(bulkModulus, massDensity, specificThermalEnergy);
}

/* ACCESSORS */
template<typename Dimension>
bool
HelmholtzEquationOfState<Dimension>::
getUpdateStatus() const {
  return needUpdate;
}

template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setUpdateStatus(bool bSet){
  needUpdate = bSet;
}
  
//------------------------------------------------------------------------------
// Set entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
      
  storeFields(massDensity, specificThermalEnergy);
      
  /* what follows is a method of breaking up the input to the fortran routine into smaller chunks */
      
  int npart = massDensity.numElements();
  int nblock = NBLOCK;
  int nrest = npart % nblock;
  int nloop = (npart - nrest)/nblock;
  unsigned int k = 0;
      
  if(needUpdate){
    for (int j=0; j<nloop; ++j)
    {
      k = j*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nblock, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
    /* now do the rest */
    if (nrest > 0)
    {
      k = nloop*nblock;
      Fortran2(wrapper_invert_helm_ed)(&nrest, &(myMassDensity->at(k)), &(mySpecificThermalEnergy->at(k)),
                                       &(myAbar->at(k)), &(myZbar->at(k)), &(myTemperature->at(k)),
                                       &(myPressure->at(k)), &mTmin, &(mySoundSpeed->at(k)), &(myGamma->at(k)),
                                       &(myEntropy->at(k)));
    }
          
  }
      
  for (int i = 0; i != npart; ++i) {
    entropy(i) = (*myEntropy)[i] / mEnergyinergpg;
  }
}


//------------------------------------------------------------------------------
// Access abar
//------------------------------------------------------------------------------
template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
HelmholtzEquationOfState<Dimension>::
abar() const {
  //return mabar;
  return *myAbar;
}

//------------------------------------------------------------------------------
// Access zbar
//------------------------------------------------------------------------------
template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
HelmholtzEquationOfState<Dimension>::
zbar() const {
  //return mzbar;
  return *myZbar;
}
      
//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
HelmholtzEquationOfState<Dimension>::valid() const {
  return (1.0);
}
  
//------------------------------------------------------------------------------
// Store Fields to local memory
//------------------------------------------------------------------------------
template<typename Dimension>
void
HelmholtzEquationOfState<Dimension>::
storeFields(const Field<Dimension, Scalar>& thisMassDensity, const Field<Dimension, Scalar>& thisSpecificThermalEnergy) const
{
      
  myMassDensity           = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>(thisMassDensity));
  mySpecificThermalEnergy = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>(thisSpecificThermalEnergy));
  myTemperature           = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>("helmTemperature",thisMassDensity.nodeList(),1.0e3*mTmin));
  myPressure              = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>("helmPressure",thisMassDensity.nodeList()));
  mySoundSpeed            = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>("helmSoundSpeed",thisMassDensity.nodeList()));
  myGamma                 = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>("helmGamma",thisMassDensity.nodeList()));
  myAbar                  = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>("helmAbar",thisMassDensity.nodeList(),mabar0));
  myZbar                  = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>("helmZbar",thisMassDensity.nodeList(),mzbar0));
  myEntropy               = shared_ptr<Field<Dimension, Scalar> >(new Field<Dimension, Scalar>("helmEntropy",thisMassDensity.nodeList()));
      
  for(unsigned int i=0; i!=myMassDensity->numElements(); ++i)
  {
    (*myMassDensity)[i]             = (*myMassDensity)[i] * mDensingpccm;
    (*mySpecificThermalEnergy)[i]   = (*mySpecificThermalEnergy)[i] * mEnergyinergpg;
          
    (*myMassDensity)[i] = ((*myMassDensity)[i]>1.0e-10 ? (*myMassDensity)[i] : 1.0e-10);
  }
}

}

