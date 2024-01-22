#ifndef MHD_CONDUCTINGFLUIDNODELIST_HH
#define MHD_CONDUCTINGFLUIDNODELIST_HH

#include "NodeList/SphNodeList.hh"
#include "Geometry/Dimension.hh"

#include <memory>

namespace Spheral {

// Predeclarations.

template <typename Dimension, typename Value> class UpdatePolicyBase;

template <typename Dimension, typename Value> class Field;
template <typename Dimension> class TableKernel;
class FileIO;

class ConductingFluidNodeList: public SphNodeList<Dim<3> >
{

   public:

   //--------------------------- Public Interface ---------------------------//
   typedef Dim<3>::Scalar Scalar;
   typedef Dim<3>::Vector Vector;
   typedef Dim<3>::Tensor Tensor;
   typedef Dim<3>::SymTensor SymTensor;

   typedef SphNodeList<Dim<3> >::FieldBaseIterator FieldBaseIterator;
   typedef SphNodeList<Dim<3> >::const_FieldBaseIterator const_FieldBaseIterator;

   //! Construct a conducting fluid node list with a constant, unchanging
   //! electrical resistivity.
   //! \param name The name by which this node list will be identified.
   //! \param eos The equation of state for the fluid.
   //! \param W The kernel function for the fluid.
   //! \param WQ The kernel function for the artificial viscosity.
   //! \param numInternalNodes The number of internal nodes in this list.
   //! \param numGhostNodes The number of ghost nodes in this list.
   ConductingFluidNodeList(const std::string& name,
      EquationOfState<Dim<3> >& eos,
      TableKernel<Dim<3> >& W,
      TableKernel<Dim<3> >& WQ,
      int numInternalNodes = 0,
      int numGhostNodes = 0);

   //! Construct a fluid node list with Spitzer resistivity and no 
   //! maximum resistivity.
   //! \param name The name by which this node list will be identified.
   //! \param eos The equation of state for the fluid.
   //! \param W The kernel function for the fluid.
   //! \param WQ The kernel function for the artificial viscosity.
   //! \param C The Spitzer model constant.
   //! \param numInternalNodes The number of internal nodes in this list.
   //! \param numGhostNodes The number of ghost nodes in this list.
   ConductingFluidNodeList(const std::string& name,
      EquationOfState<Dim<3> >& eos,
      TableKernel<Dim<3> >& W,
      TableKernel<Dim<3> >& WQ,
      const Scalar& C,
      int numInternalNodes = 0,
      int numGhostNodes = 0);

   //! Construct a fluid node list with Spitzer resistivity.
   //! \param name The name by which this node list will be identified.
   //! \param eos The equation of state for the fluid.
   //! \param W The kernel function for the fluid.
   //! \param WQ The kernel function for the artificial viscosity.
   //! \param C The Spitzer model constant.
   //! \param Rmax The maximum parallel resistivity delivered by the Spitzer model.
   //! \param numInternalNodes The number of internal nodes in this list.
   //! \param numGhost TheNodes number of ghost nodes in this list.
   ConductingFluidNodeList(const std::string& name,
      EquationOfState<Dim<3> >& eos,
      TableKernel<Dim<3> >& W,
      TableKernel<Dim<3> >& WQ,
      double C,
      double Rmax,
      int numInternalNodes = 0,
      int numGhostNodes = 0);

   //! Construct a fluid node list with the modified Spitzer resistivity, 
   //! accounting for Anomolous effects at the given minimum density threshhold.
   //! \param name The name by which this node list will be identified.
   //! \param eos The equation of state for the fluid.
   //! \param W The kernel function for the fluid.
   //! \param WQ The kernel function for the artificial viscosity.
   //! \param C The Spitzer model constant.
   //! \param Rmax The maximum parallel resistivity delivered by the Spitzer model.
   //! \param rhoMin The minimum density threshhold for charge carriers.
   //! \param numInternalNodes The number of internal nodes in this list.
   //! \param numGhostNodes The number of ghost nodes in this list.
   ConductingFluidNodeList(const std::string& name,
      EquationOfState<Dim<3> >& eos,
      TableKernel<Dim<3> >& W,
      TableKernel<Dim<3> >& WQ,
      double C,
      double Rmax,
      double rhoMin,
      int numInternalNodes = 0,
      int numGhostNodes = 0);

   //! Destructor
   ~ConductingFluidNodeList();

   //! Access the fluid's (non-const) magnetic induction.
   Field<Dim<3>, Vector>& magneticInduction();

   //! Access the fluid's (const) magnetic induction.
   const Field<Dim<3>, Vector>& magneticInduction() const;
 
   //! Access the fluid's (non-const) magnetic induction.
   Field<Dim<3>, Scalar>& magneticDivergence();

   //! Access the fluid's (const) magnetic induction.
   const Field<Dim<3>, Scalar>& magneticDivergence() const;
 
   //! Access the (non-const) time derivative of the fluid's magnetic induction.
   Field<Dim<3>, Vector>& DBDt();
   
   //! Access the (const) time derivative of the fluid's magnetic induction.
   const Field<Dim<3>, Vector>& DBDt() const;
   
   //! Access the fluid's (non-const) current density.
   Field<Dim<3>, Vector>& currentDensity();

   //! Access the fluid's (const) current density.
   const Field<Dim<3>, Vector>& currentDensity() const;
   
   //! Access the fluid's (non-const) resistivity.
   Field<Dim<3>, Scalar>& resistivity();

   //! Access the fluid's (const) resistivity.
   const Field<Dim<3>, Scalar>& resistivity() const;

   //! Access the fluid's (non-const) total specific energy.
   Field<Dim<3>, Scalar>& totalSpecificEnergy();

   //! Access the fluid's (const) total specific energy.
   const Field<Dim<3>, Scalar>& totalSpecificEnergy() const;
 
   //! Access the (non-const) time derivative of the fluid's total specific energy.
   Field<Dim<3>, Scalar>& DeDt();
   
   //! Access the (const) time derivative of the fluid's total specific energy.
   const Field<Dim<3>, Scalar>& DeDt() const;
   
   //! Access the (non-const) time derivative of the fluid's magnetic induction.
   //! Access the fluid's (const) resistivity update policy.
   const std::shared_ptr<UpdatePolicyBase<Dim<3>, Field<Dim<3>, Scalar> > >& resistivityPolicy() const;

   //! Access the fluid's (non-const) resistivity update policy.
   std::shared_ptr<UpdatePolicyBase<Dim<3>, Field<Dim<3>, Scalar> > >& resistivityPolicy();

   //! Dump the state to a restart file.
   void dumpState(FileIO& file, const std::string& pathName) const;

   //! Restore the state from a restart file.
   void restoreState(const FileIO& file, const std::string& pathName);

   private:

   // No copy constructor or assignment operator.
   ConductingFluidNodeList(const ConductingFluidNodeList&);
   ConductingFluidNodeList& operator=(const ConductingFluidNodeList&);

   // Magnetic induction, time derivative, current density.
   std::shared_ptr<Field<Dim<3>, Vector> > 
      mMagneticInduction, mDBDt, mCurrentDensity;

   // Resistivity and energy fields.
   std::shared_ptr<Field<Dim<3> , Scalar> > mDivB, mResistivity, 
                                              mTotalSpecificEnergy, mDeDt;

   // Resistivity update policy.
   std::shared_ptr<UpdatePolicyBase<Dim<3>, Field<Dim<3>, Scalar> > > mResistivityPolicy;

}; // end class ConductingFluidNodeList

} // end namespace Spheral

#endif
