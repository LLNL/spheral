#include "NodeList/SphNodeList.hh"
#include "MHD/ConductingFluidNodeList.hh"
#include "Kernel/TableKernel.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/NonDynamicState.hh"
#include "MHD/SpitzerResistivityUpdatePolicy.hh"
#include "FileIO/FileIO.hh"
#include "MHD/MHDFieldNames.hh"

namespace Spheral {

//----------------------------------------------------------------------------
ConductingFluidNodeList::
ConductingFluidNodeList(const std::string& name,
                        EquationOfState<Dim<3> >& eos,
                        TableKernel<Dim<3> >& W,
                        TableKernel<Dim<3> >& WQ,
                        int numInternalNodes,
                        int numGhostNodes):
   SphNodeList<Dim<3> >::SphNodeList(name, eos, W, WQ, numInternalNodes, numGhostNodes),
   mMagneticInduction(new Field<Dim<3>, Vector>(MHDFieldNames::magneticInduction, *this)),
   mDBDt(new Field<Dim<3>, Vector>(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + MHDFieldNames::magneticInduction, *this)),
   mCurrentDensity(new Field<Dim<3>, Vector>(MHDFieldNames::currentDensity, *this)),
   mDivB(new Field<Dim<3>, Scalar>(MHDFieldNames::magneticDivergence, *this)),
        mResistivity(new Field<Dim<3>, Scalar>(MHDFieldNames::resistivity, *this)),
        mTotalSpecificEnergy(new Field<Dim<3>, Scalar>(MHDFieldNames::totalSpecificEnergy, *this)),
        mDeDt(new Field<Dim<3>, Dim<3>::Scalar>(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + MHDFieldNames::totalSpecificEnergy, *this)),
   mResistivityPolicy(new NonDynamicState<Dim<3>, Field<Dim<3>, Scalar> >())
{
} // end constructor
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
ConductingFluidNodeList::
ConductingFluidNodeList(const std::string& name,
                        EquationOfState<Dim<3> >& eos,
                        TableKernel<Dim<3> >& W,
                        TableKernel<Dim<3> >& WQ,
                        const Scalar& C,
                        int numInternalNodes,
                        int numGhostNodes):
   SphNodeList<Dim<3> >::SphNodeList(name, eos, W, WQ, numInternalNodes, numGhostNodes),
   mMagneticInduction(new Field<Dim<3>, Vector>(MHDFieldNames::magneticInduction, *this)),
   mDBDt(new Field<Dim<3>, Vector>(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + MHDFieldNames::magneticInduction, *this)),
   mCurrentDensity(new Field<Dim<3>, Vector>(MHDFieldNames::currentDensity, *this)),
   mDivB(new Field<Dim<3>, Scalar>(MHDFieldNames::magneticDivergence, *this)),
        mResistivity(new Field<Dim<3>, Dim<3>::Scalar>(MHDFieldNames::resistivity, *this)),
        mTotalSpecificEnergy(new Field<Dim<3>, Dim<3>::Scalar>(MHDFieldNames::totalSpecificEnergy, *this)),
        mDeDt(new Field<Dim<3>, Dim<3>::Scalar>(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + MHDFieldNames::totalSpecificEnergy, *this)),
   mResistivityPolicy(new SpitzerResistivityUpdatePolicy(C))
{
} // end constructor
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
ConductingFluidNodeList::
ConductingFluidNodeList(const std::string& name,
                        EquationOfState<Dim<3> >& eos,
                        TableKernel<Dim<3> >& W,
                        TableKernel<Dim<3> >& WQ,
                        double C,
                        double Rmax,
                        int numInternalNodes,
                        int numGhostNodes):
   SphNodeList<Dim<3> >::SphNodeList(name, eos, W, WQ, numInternalNodes, numGhostNodes),
   mMagneticInduction(new Field<Dim<3>, Vector>(MHDFieldNames::magneticInduction, *this)),
   mDBDt(new Field<Dim<3>, Vector>(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + MHDFieldNames::magneticInduction, *this)),
   mCurrentDensity(new Field<Dim<3>, Vector>(MHDFieldNames::currentDensity, *this)),
   mDivB(new Field<Dim<3>, Scalar>(MHDFieldNames::magneticDivergence, *this)),
        mResistivity(new Field<Dim<3>, Dim<3>::Scalar>(MHDFieldNames::resistivity, *this)),
        mTotalSpecificEnergy(new Field<Dim<3>, Dim<3>::Scalar>(MHDFieldNames::totalSpecificEnergy, *this)),
        mDeDt(new Field<Dim<3>, Dim<3>::Scalar>(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + MHDFieldNames::totalSpecificEnergy, *this)),
   mResistivityPolicy(new SpitzerResistivityUpdatePolicy(C, Rmax))
{
} // end constructor
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
ConductingFluidNodeList::
ConductingFluidNodeList(const std::string& name,
                        EquationOfState<Dim<3> >& eos,
                        TableKernel<Dim<3> >& W,
                        TableKernel<Dim<3> >& WQ,
                        double C,
                        double Rmax,
                        double rhoMin,
                        int numInternalNodes,
                        int numGhostNodes):
   SphNodeList<Dim<3> >::SphNodeList(name, eos, W, WQ, numInternalNodes, numGhostNodes),
   mMagneticInduction(new Field<Dim<3>, Vector>(MHDFieldNames::magneticInduction, *this)),
   mDBDt(new Field<Dim<3>, Vector>(IncrementState<Dim<3>, Field<Dim<3>, Vector> >::prefix() + MHDFieldNames::magneticInduction, *this)),
   mCurrentDensity(new Field<Dim<3>, Vector>(MHDFieldNames::currentDensity, *this)),
   mDivB(new Field<Dim<3>, Scalar>(MHDFieldNames::magneticDivergence, *this)),
        mResistivity(new Field<Dim<3>, Dim<3>::Scalar>(MHDFieldNames::resistivity, *this)),
        mTotalSpecificEnergy(new Field<Dim<3>, Dim<3>::Scalar>(MHDFieldNames::totalSpecificEnergy, *this)),
        mDeDt(new Field<Dim<3>, Dim<3>::Scalar>(IncrementState<Dim<3>, Field<Dim<3>, Scalar> >::prefix() + MHDFieldNames::totalSpecificEnergy, *this)),
   mResistivityPolicy(new SpitzerResistivityUpdatePolicy(C, Rmax, rhoMin))
{
} // end constructor
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
ConductingFluidNodeList::
~ConductingFluidNodeList()
{
} // end destructor
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
const Field<Dim<3>, Dim<3>::Vector>& 
ConductingFluidNodeList::
magneticInduction() const
{
   return *mMagneticInduction;
} // end magneticInduction
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
Field<Dim<3>, Dim<3>::Vector>& 
ConductingFluidNodeList::
magneticInduction()
{
   return *mMagneticInduction;
} // end magneticInduction
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
const Field<Dim<3>, Dim<3>::Scalar>& 
ConductingFluidNodeList::
magneticDivergence() const
{
   return *mDivB;
} // end magneticDivergence
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
Field<Dim<3>, Dim<3>::Scalar>& 
ConductingFluidNodeList::
magneticDivergence()
{
   return *mDivB;
} // end magneticDivergence
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
const Field<Dim<3>, Dim<3>::Vector>& 
ConductingFluidNodeList::
DBDt() const
{
   return *mDBDt;
} // end DBDt
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
Field<Dim<3>, Dim<3>::Vector>& 
ConductingFluidNodeList::
DBDt()
{
   return *mDBDt;
} // end DBDt
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
const Field<Dim<3>, Dim<3>::Vector>& 
ConductingFluidNodeList::
currentDensity() const
{
   return *mCurrentDensity;
} // end currentDensity
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
Field<Dim<3>, Dim<3>::Vector>& 
ConductingFluidNodeList::
currentDensity()
{
   return *mCurrentDensity;
} // end currentDensity
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
const Field<Dim<3>, Dim<3>::Scalar>& 
ConductingFluidNodeList::
totalSpecificEnergy() const
{
   return *mTotalSpecificEnergy;
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
Field<Dim<3>, Dim<3>::Scalar>& 
ConductingFluidNodeList::
totalSpecificEnergy()
{
   return *mTotalSpecificEnergy;
} 
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
const Field<Dim<3>, Dim<3>::Scalar>& 
ConductingFluidNodeList::
DeDt() const
{
   return *mDeDt;
} // end currentDensity
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
Field<Dim<3>, Dim<3>::Scalar>& 
ConductingFluidNodeList::
DeDt()
{
   return *mDeDt;
} // end currentDensity
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
Field<Dim<3>, Dim<3>::Scalar>& 
ConductingFluidNodeList::
resistivity()
{
   return *mResistivity;
} // end resistivity
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
const Field<Dim<3>, Dim<3>::Scalar>& 
ConductingFluidNodeList::
resistivity() const
{
   return *mResistivity;
} // end resistivity
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
const std::shared_ptr<UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> > >&
ConductingFluidNodeList::
resistivityPolicy() const
{
   return mResistivityPolicy;
} // end resistivityPolicy
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
std::shared_ptr<UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> > >&
ConductingFluidNodeList::
resistivityPolicy() 
{
   return mResistivityPolicy;
} // end resistivityPolicy
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void
ConductingFluidNodeList::
dumpState(FileIO& file, const std::string& pathName) const {

  // Dump the ancestor class.
  SphNodeList<Dim<3> >::dumpState(file, pathName);

  // Dump each of the internal fields of the ConductingFluidNodeList.
  file.write(*mMagneticInduction, pathName + "/magneticInduction");
  file.write(*mDBDt, pathName + "/DBDt");
  file.write(*mCurrentDensity, pathName + "/currentDensity");
  file.write(*mDivB, pathName + "/divB");
  file.write(*mResistivity, pathName + "/resistivity");
  file.write(*mTotalSpecificEnergy, pathName + "/totalSpecificEnergy");
  file.write(*mDeDt, pathName + "/DeDt");

  // Now etch our resistivity policy in there.
  int RPolicyType;
  if (dynamic_cast<NonDynamicState<Dim<3>, Field<Dim<3>, Scalar> >* const>(mResistivityPolicy.get()) != 0)
  {
    RPolicyType = 0;
    file.write(RPolicyType, pathName + "/resistivityPolicy");
  }
  else // Spitzer policy
  {
    SpitzerResistivityUpdatePolicy* const policy = 
      dynamic_cast<SpitzerResistivityUpdatePolicy* const>(mResistivityPolicy.get());
    CHECK(policy != 0);
    RPolicyType = 1;
    file.write(RPolicyType, pathName + "/resistivityPolicy");
    file.write(policy->C(), pathName + "/SpitzerResistivity_C");
    file.write(policy->Rmax(), pathName + "/SpitzerResistivity_Rmax");
    file.write(policy->rhoMin(), pathName + "/SpitzerResistivity_rhoMin");
  }
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
void
ConductingFluidNodeList::
restoreState(const FileIO& file, const std::string& pathName) {

  // Restore the ancestor class.
  SphNodeList<Dim<3> >::restoreState(file, pathName);

  // Restore each of the internal fields of the FluidNodeList.
  file.read(*mMagneticInduction, pathName + "/magneticInduction");
  file.read(*mDBDt, pathName + "/DBDt");
  file.read(*mDivB, pathName + "/divB");
  file.read(*mCurrentDensity, pathName + "/currentDensity");
  file.read(*mResistivity, pathName + "/resistivity");
  file.read(*mTotalSpecificEnergy, pathName + "/totalSpecificEnergy");
  file.read(*mDeDt, pathName + "/DeDt");

  // Resistivity policy.
  int RpolicyType;
  file.read(RpolicyType, pathName + "/resistivityPolicy");

  if (RpolicyType == 0) // Constant
  {
    mResistivityPolicy.reset(new NonDynamicState<Dim<3>, Field<Dim<3>, Scalar> >());
  }
  else // Spitzer policy
  {
    CHECK(RpolicyType == 1) // Spitzer
    double C, Rmax, rhoMin;
    file.read(C, pathName + "/SpitzerResistivity_C");
    file.read(Rmax, pathName + "/SpitzerResistivity_Rmax");
    file.read(rhoMin, pathName + "/SpitzerResistivity_rhoMin");
    mResistivityPolicy.reset(new SpitzerResistivityUpdatePolicy(C, Rmax, rhoMin));
  }
}  

//----------------------------------------------------------------------------

} // end namespace Spheral
