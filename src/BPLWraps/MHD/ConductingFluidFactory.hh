#include "MHD/ConductingFluidNodeList.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {
namespace MHDSpace {


class ConductingFluidFactory
{
   public:

   ConductingFluidFactory() {}
   ~ConductingFluidFactory() {}

   //! Here's a ConductingFluidNodeList constructor for a perfectly-conducting 
   //! fluid.
   boost::shared_ptr<ConductingFluidNodeList>
   PerfectlyConductingFluid(const std::string& name, 
                            EquationOfState<Dim<3> >& eos, 
                            KernelSpace::TableKernel<Dim<3> >& W,
                            KernelSpace::TableKernel<Dim<3> >& WQ,
                            double Z, 
                            int numInternalNodes = 0, 
                            int numGhostNodes = 0)
   {
      boost::shared_ptr<ConductingFluidNodeList>
         fluid(new ConductingFluidNodeList(name, eos, W, WQ, Z, numInternalNodes, numGhostNodes));
      return fluid;
   } // end PerfectlyConductingFluid

   //! Here's a ConductingFluidNodeList with constant scalar conductivity.
   boost::shared_ptr<ConductingFluidNodeList<Dimension> >
   ResistiveFluid(const std::string& name, 
                  EquationOfState<Dimension>& eos, 
                  KernelSpace::TableKernel<Dimension>& W,
                  KernelSpace::TableKernel<Dimension>& WQ,
                  double R0,
                  double Z, 
                  int numInternalNodes = 0, 
                  int numGhostNodes = 0)
   {
      boost::shared_ptr<ConductingFluidNodeList<Dimension> >
         fluid(new ConductingFluidNodeList<Dimension>(name, eos, W, WQ, R0, Z, numInternalNodes, numGhostNodes));
      return fluid;
   } // end ResistiveFluid


   //! Here's a ConductingFluidNodeList with constant tensor conductivity.
   boost::shared_ptr<ConductingFluidNodeList<Dimension> >
   ResistiveFluid(const std::string& name, 
                  EquationOfState<Dimension>& eos, 
                  KernelSpace::TableKernel<Dimension>& W,
                  KernelSpace::TableKernel<Dimension>& WQ,
                  const typename Dimension::SymTensor& R0,
                  double Z, 
                  int numInternalNodes = 0, 
                  int numGhostNodes = 0)
   {
      boost::shared_ptr<ConductingFluidNodeList<Dimension> >
         fluid(new ConductingFluidNodeList<Dimension>(name, eos, W, WQ, R0, Z, numInternalNodes, numGhostNodes));
      return fluid;
   } // end ResistiveFluid


   //! Here's a ConductingFluidNodeList with Spitzer conductivity.
   boost::shared_ptr<ConductingFluidNodeList<Dimension> >
   SpitzerResistiveFluid(const std::string& name, 
                         EquationOfState<Dimension>& eos, 
                         KernelSpace::TableKernel<Dimension>& W,
                         KernelSpace::TableKernel<Dimension>& WQ,
                         double C,
                         double Rmax,
                         double Z, 
                         int numInternalNodes = 0, 
                         int numGhostNodes = 0)
   {
      boost::shared_ptr<ConductingFluidNodeList<Dimension> >
         fluid(new ConductingFluidNodeList<Dimension>(name, eos, W, WQ, C, Rmax, Z, numInternalNodes, numGhostNodes));
      return fluid;
   } // end SpitzerResistiveFluid


   //! Here's a ConductingFluidNodeList with Density-limited Spitzer conductivity.
   boost::shared_ptr<ConductingFluidNodeList<Dimension> >
   DensityLimitedSpitzerResistiveFluid(const std::string& name, 
                                       EquationOfState<Dimension>& eos, 
                                       KernelSpace::TableKernel<Dimension>& W,
                                       KernelSpace::TableKernel<Dimension>& WQ,
                                       double C,
                                       double Rmax,
                                       double rhoMin,
                                       double Z, 
                                       int numInternalNodes = 0, 
                                       int numGhostNodes = 0)
   {
      boost::shared_ptr<ConductingFluidNodeList<Dimension> >
         fluid(new ConductingFluidNodeList<Dimension>(name, eos, W, WQ, C, Rmax, rhoMin, Z, numInternalNodes, numGhostNodes));
      return fluid;
   } // end SpitzerResistiveFluid

   private:

   ConductingFluidFactory(const ConductingFluidFactory&);
   ConductingFluidFactory& operator=(const ConductingFluidFactory&);

}; // end class ConductingFluidFactory


}
}
