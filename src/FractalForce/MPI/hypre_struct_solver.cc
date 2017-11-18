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
  void hypre_struct_solver(vector <Point*>& p_points_left,vector <Point*>& p_points_right,Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    const unsigned int zoom=Misc::pow(2,fractal.get_level_max()-level);
    const unsigned int ndim=3;
    const unsigned int stencil_size=7;
    int length=fractal.get_grid_length();
    vector <int>pl(ndim);
    vector <int>pr(ndim);
    vector <int>dprl(ndim);
    int apl[3];
    int apr[3];
    vector <int>bplr(ndim);
    int offsets[7][3]={{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    HYPRE_StructGrid* p_hypre_grid;
    HYPRE_StructGrid& hypre_grid=*p_hypre_grid;
    HYPRE_StructGridCreate(mem.p_mess->FractalWorld,ndim,p_hypre_grid); 
    unsigned int largest=-1;
    for(unsigned int np=0;np<p_points_left.size();np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	unsigned int how_large=1;
	for(unsigned int ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    how_large*=apr[ni]-apl[ni]+1;
	  }
	largest=max(largest,how_large);
	HYPRE_StructGridSetExtents(hypre_grid,apl,apr);
      }
    int periods[3];
    if(fractal.get_periodic())
      {
	periods[0]=0;
	periods[1]=0;
	periods[2]=0;
      }
    else
      {
	periods[0]=length;
	periods[1]=length;
	periods[2]=length;
      }
    HYPRE_StructGridSetPeriodic(hypre_grid,periods);
    HYPRE_StructGridAssemble(hypre_grid);
    //
    HYPRE_StructStencil* p_hypre_stencil;
    HYPRE_StructStencil& hypre_stencil=*p_hypre_stencil;
    HYPRE_StructStencilCreate(ndim,stencil_size,p_hypre_stencil);
    for(unsigned int stencil_entry=0;stencil_entry<stencil_size;stencil_entry++)
      HYPRE_StructStencilSetElement(hypre_stencil,stencil_entry,offsets[stencil_entry]);
    //
    HYPRE_StructMatrix* p_hypre_matrix;
    HYPRE_StructMatrix& hypre_matrix= *p_hypre_matrix;
    HYPRE_StructMatrixCreate(mem.p_mess->FractalWorld,hypre_grid,hypre_stencil,p_hypre_matrix);
    HYPRE_StructMatrixInitialize(hypre_matrix);
    int stencil_indices[7]={0,1,2,3,4,5,6};
    vector <double>values(largest*7);
    unsigned n=0;
    for(unsigned int p=0;p<largest;p++)
      {
	values[n]=-6.0;
	n++;
	for(int q=0;q<6;q++)
	  {
	    values[n]=1.0;
	    n++;
	  }
      }
    for(unsigned int np=0;np<p_points_left.size();np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	unsigned int how_large=1;
	for(unsigned int ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    how_large*=apr[ni]-apl[ni]+1;
	  }
	HYPRE_StructMatrixSetBoxValues(hypre_matrix, apl, apr, how_large*7,stencil_indices, values);
      }
    //
    // something for the boundaries
    values.assign(7,0.0);
    vector < vector <bool> > bcs(largest);
    for(unsigned int np=0;np<p_points_left.size();np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	for(unsigned int ni=0;ni<ndim;ni++)
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
			HYPRE_StructMatrixSetBoxValues(hypre_matrix, bplr, bplr, outs,stencil_indices, values);
		      }
		    ni++;
		  }
	      }
	  }
      }
    bcs.clear();
    HYPRE_StructMatrixAssemble(hypre_matrix);
    //
    HYPRE_StructVector* p_pots;
    HYPRE_StructVector& pots=*p_pots;
    HYPRE_StructVectorCreate(mem.p_mess->FractalWorld,hypre_grid,p_pots);
    HYPRE_StructVectorInitialize(pots);
    HYPRE_StructVector* p_dens;
    HYPRE_StructVector& dens=*p_dens;
    HYPRE_StructVectorCreate(mem.p_mess->FractalWorld,hypre_grid,p_dens);
    HYPRE_StructVectorInitialize(dens);
    //
    const double pi = 4.0*atan(1.0);
    double g_c=4.0*pi/static_cast<double>(length*length)*pow(4.0,-level);
    //
    for(unsigned int np=0;np<p_points_left.size();np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	for(unsigned int ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    dprl[ni]=apr[ni]-apl[ni]+1;
	  }
	vector <double>potss(largest);
	vector <double>denss(largest);
	p_points_left[np]->get_potss_denss(dprl,g_c,potss,denss);
	HYPRE_StructVectorSetBoxValues(pots,apl,apr,potss);
	HYPRE_StructVectorSetBoxValues(dens,apl,apr,denss);
      }
    HYPRE_StructVectorAssemble(pots);
    HYPRE_StructVectorAssemble(dens);

    HYPRE_StructSolver* p_hypre_solver;
    HYPRE_StructSolver& hypre_solver= *p_hypre_solver;

    if(mem.hypre_solver == "PCG")
      {
	HYPRE_StructPCGCreate(mem.p_mess->FractalWorld,p_hypre_solver);
	HYPRE_StructPCGSetTol(hypre_solver,fractal.get_epsilon_sor());
	HYPRE_StructPCGSetMaxIter(hypre_solver,fractal.get_maxits());
	HYPRE_StructPCGSetPrintLevel(hypre_solver,2);
	HYPRE_StructPCGSetTwoNorm(hypre_solver, 1 );
	HYPRE_StructPCGSetRelChange(hypre_solver, 0 );
	HYPRE_StructPCGSetLogging(hypre_solver, 1);
      }
    else
      {
	cout << "solver set " << mem.hypre_solver << endl;
	assert(false);
      }
    if(mem.hypre_precond == "SMG")
      {
	/* Use symmetric SMG as preconditioner */
	HYPRE_StructSolver* p_hypre_precond;
	HYPRE_StructSolver& hypre_precond= *p_hypre_precond;
	HYPRE_StructSMGCreate(mem.p_mess->FractalWorld, p_hypre_precond);
	HYPRE_StructSMGSetMemoryUse(hypre_precond, 0);
	HYPRE_StructSMGSetMaxIter(hypre_precond, 1);
	HYPRE_StructSMGSetTol(hypre_precond, 0.0);
	HYPRE_StructSMGSetNumPreRelax(hypre_precond, 1);
	HYPRE_StructSMGSetNumPostRelax(hypre_precond, 1);
      /* Set the preconditioner and solve */
	if(mem.hypre_solver == "PCG")
	  HYPRE_StructPCGSetPrecond(hypre_solver, HYPRE_StructSMGSolve,
				    HYPRE_StructSMGSetup, hypre_precond);
      }
    if(mem.hypre_solver == "PCG")
      {
	HYPRE_StructPCGSetup(hypre_solver,hypre_matrix,dens,pots);
	HYPRE_StructPCGSolve(hypre_solver,hypre_matrix,dens,pots);
	HYPRE_StructPCGDestroy(hypre_precond);
	HYPRE_StructPCGDestroy(hypre_solver);
      }
    HYPRE_StructGridDestroy(hypre_grid);
    HYPRE_StructStencilDestroy(hypre_stencil);
    HYPRE_StructMatrixDestroy(hypre_matrix);
    HYPRE_StructVectorDestroy(dens);
    //
    for(unsigned int np=0;np<p_points_left.size();np++)
      {
	p_points_left[np]->get_pos_point(pl);
	p_points_right[np]->get_pos_point(pr);
	unsigned int how_large=1;
	for(unsigned int ni=0;ni<ndim;ni++)
	  {
	    apl[ni]=pl[ni]/zoom;
	    apr[ni]=pr[ni]/zoom;
	    dprl[ni]=apr[ni]-apl[ni]+1;
	    how_large*=dprl[ni];
	  }
	HYPRE_StructVectorGetBoxValues(pots,apl,apr,potss);
	p_points_left[np]->set_potss(dprl,potss);
      }
    //
    HYPRE_StructVectorDestroy(pots);
  }
}
