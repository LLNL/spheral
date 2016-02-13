#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#ifndef _Hypre_Defined_
#define _Hypre_Defined_
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
#endif
namespace FractalSpace
{
  void hypre_struct_solver(vector <Point*>& p_points_left,vector <Point*>& p_points_right,Fractal& fractal,Fractal_Memory& mem,
			   const int& level,const bool& buffer_groups)
  {
    typedef unsigned int uint;
    MPI_Comm HypreComm;
    if(buffer_groups)
      HypreComm=mem.p_mess->HypreWorld;
    else
      HypreComm=mem.p_mess->SELF;
    ofstream& FileHypre=mem.p_file->FileHypre;
    const uint zoom=Misc::pow(2,fractal.get_level_max()-level);
    const uint ndim=3;
    const uint stencil_size=7;
    const uint total_boxes=p_points_left.size();
    const uint length=fractal.get_grid_length();
    const uint length_box=length*Misc::pow(2,level);
    vector <int>pl(ndim);
    vector <int>pr(ndim);
    int dprl[3];
    int apl[3];
    int apr[3];
    int bplr[3];
    int offsets[7][3]={{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    HYPRE_StructGrid hypre_grid;
    HYPRE_StructGridCreate(HypreComm,ndim,&hypre_grid); 
    uint largest=-1;
    for(uint np=0;np<total_boxes;np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	uint how_large=1;
	for(uint ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    how_large*=apr[ni]-apl[ni]+1;
	  }
	largest=max(largest,how_large);
	HYPRE_StructGridSetExtents(hypre_grid,apl,apr);
      }
    int periods[3]={0,0,0};
    if(fractal.get_periodic())
      {
	periods[0]=length_box;
	periods[1]=length_box;
	periods[2]=length_box;
      }
    HYPRE_StructGridSetPeriodic(hypre_grid,periods);
    HYPRE_StructGridAssemble(hypre_grid);
    //
    HYPRE_StructStencil hypre_stencil;
    HYPRE_StructStencilCreate(ndim,stencil_size,&hypre_stencil);
    for(uint stencil_entry=0;stencil_entry<stencil_size;stencil_entry++)
      HYPRE_StructStencilSetElement(hypre_stencil,stencil_entry,offsets[stencil_entry]);
    //
    HYPRE_StructMatrix hypre_matrix;
    HYPRE_StructMatrixCreate(HypreComm,hypre_grid,hypre_stencil,&hypre_matrix);
    HYPRE_StructMatrixInitialize(hypre_matrix);
    int stencil_indices[7]={0,1,2,3,4,5,6};
    double* Mvalues=new double[largest*7];
    uint n=0;
    for(uint p=0;p<largest;p++)
      {
	Mvalues[n]=-6.0;
	n++;
	for(int q=0;q<6;q++)
	  {
	    Mvalues[n]=1.0;
	    n++;
	  }
      }
    for(uint np=0;np<total_boxes;np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	uint how_large=1;
	for(uint ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    how_large*=apr[ni]-apl[ni]+1;
	  }
	HYPRE_StructMatrixSetBoxValues(hypre_matrix, apl, apr, how_large*7,stencil_indices, Mvalues);
      }
    //
    // something for the boundaries
    for(int ni=0;ni<7;ni++)
      Mvalues[ni]=0;
    vector < vector <bool> > bcs(largest);
    for(uint np=0;np<total_boxes;np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	for(uint ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    dprl[ni]=apr[ni]-apl[ni]+1;
	  }
	if(!p_points_left[np]->get_bcs(bcs,dprl))
	  continue;
	int ni=0;
	for(int ni2=apl[2];ni2<=apr[2];ni2++)
	  {
	    bplr[2]=ni2;
	    for(int ni1=apl[1];ni1<=apr[1];ni1++)
	      {
		bplr[1]=ni1;
		for(int ni0=apl[0];ni0<=apr[0];ni0++)
		  {
		    if(bcs[ni][0])
		      {
			int outs=0;
			bplr[0]=ni0;
			for(int q=1;q<7;q++)
			  {
			    if(!bcs[ni][q])
			      continue;
			    stencil_indices[outs]=q;
			    outs++;
			  }
			HYPRE_StructMatrixSetBoxValues(hypre_matrix, bplr, bplr, outs,stencil_indices, Mvalues);
		      }
		    ni++;
		  }
	      }
	  }
      }
    delete [] Mvalues;
    bcs.clear();
    HYPRE_StructMatrixAssemble(hypre_matrix);
    //
    HYPRE_StructVector hypre_pots;
    HYPRE_StructVectorCreate(HypreComm,hypre_grid,&hypre_pots);
    HYPRE_StructVectorInitialize(hypre_pots);
    HYPRE_StructVector hypre_dens;
    HYPRE_StructVectorCreate(HypreComm,hypre_grid,&hypre_dens);
    HYPRE_StructVectorInitialize(hypre_dens);
    //
    const double pi = 4.0*atan(1.0);
    double g_c=4.0*pi/static_cast<double>(length*length)*pow(4.0,-level);
    //
    double* potss=new double[largest];
    double* denss=new double[largest];
    for(uint np=0;np<total_boxes;np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	for(uint ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    dprl[ni]=apr[ni]-apl[ni]+1;
	  }
	p_points_left[np]->get_potss_denss(dprl,g_c,potss,denss);
	HYPRE_StructVectorSetBoxValues(hypre_pots,apl,apr,potss);
	HYPRE_StructVectorSetBoxValues(hypre_dens,apl,apr,denss);
      }
    HYPRE_StructVectorAssemble(hypre_pots);
    HYPRE_StructVectorAssemble(hypre_dens);

    HYPRE_StructSolver hypre_solver;

    if(mem.hypre_solver == "PCG")
      {
	HYPRE_StructPCGCreate(HypreComm,&hypre_solver);
	HYPRE_StructPCGSetMaxIter(hypre_solver,fractal.get_maxits());
	HYPRE_StructPCGSetTol(hypre_solver,fractal.get_epsilon_sor());
	HYPRE_StructPCGSetTwoNorm(hypre_solver, 1 );
	HYPRE_StructPCGSetRelChange(hypre_solver, 0 );
	HYPRE_StructPCGSetPrintLevel(hypre_solver,2);
	HYPRE_StructPCGSetLogging(hypre_solver, 1);
      }
    else
      {
	FileHypre << "solver set " << mem.hypre_solver << endl;
	assert(false);
      }
    HYPRE_StructSolver hypre_precond;
    if(mem.hypre_precond == "SMG")
      {
	/* Use symmetric SMG as preconditioner */
	HYPRE_StructSMGCreate(HypreComm, &hypre_precond);
	HYPRE_StructSMGSetMemoryUse(hypre_precond, 0);
	HYPRE_StructSMGSetMaxIter(hypre_precond, 1);
	HYPRE_StructSMGSetTol(hypre_precond, 0.0);
	HYPRE_StructSMGSetNumPreRelax(hypre_precond, 1);
	HYPRE_StructSMGSetNumPostRelax(hypre_precond, 1);
      /* Set the preconditioner and solve */
	if(mem.hypre_solver == "PCG")
	  HYPRE_StructPCGSetPrecond(hypre_solver, HYPRE_StructSMGSolve,
				    HYPRE_StructSMGSetup,hypre_precond);
      }
    if(mem.hypre_solver == "PCG")
      {
	HYPRE_StructPCGSetup(hypre_solver,hypre_matrix,hypre_dens,hypre_pots);
	HYPRE_StructPCGSolve(hypre_solver,hypre_matrix,hypre_dens,hypre_pots);
	HYPRE_StructPCGDestroy(hypre_solver);
      }
    HYPRE_StructGridDestroy(hypre_grid);
    HYPRE_StructStencilDestroy(hypre_stencil);
    HYPRE_StructMatrixDestroy(hypre_matrix);
    HYPRE_StructVectorDestroy(hypre_dens);
    //
    for(uint np=0;np<total_boxes;np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	uint how_large=1;
	for(uint ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    dprl[ni]=apr[ni]-apl[ni]+1;
	    how_large*=dprl[ni];
	  }
	HYPRE_StructVectorGetBoxValues(hypre_pots,apl,apr,potss);
	p_points_left[np]->set_potss(dprl,potss);
      }
    //
    HYPRE_StructVectorDestroy(hypre_pots);
    delete [] potss;
    delete [] denss;
  }
}
