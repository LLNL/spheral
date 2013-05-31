#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#ifndef _Hypre_Defined_
#define _Hypre_Defined_
#include "HYPRE.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#endif
namespace FractalSpace
{
  void hypre_solver(Fractal& frac,Fractal_Memory& mem,const int& level)
  {
    ofstream& FH=mem.p_file->FileHypre;
    char describe[100];
    FH << " enter hypre solver " << level << endl;
    vector <Point*>hypre_points;
    if(!hypre_ij_numbering(mem,frac,hypre_points,level))
      {
	FH << " nothing here hypre solver " << level << endl;
	return;
      }
    FH << " really enter hypre solver a " << level << endl;
    bool inside,edge,buff,pass;
    const int FractalRank=mem.p_mess->FractalRank;
    const int HypreRank=mem.p_mess->HypreRank;
    MPI_Comm HypreComm=mem.p_mess->HypreWorld;
    HYPRE_IJMatrix ij_matrix;
    HYPRE_ParCSRMatrix par_matrix;
    const int ilower=mem.ij_offsets[HypreRank];
    const int iupper=ilower+mem.ij_counts[HypreRank]-1;
    const int jlower=ilower;
    const int jupper=iupper;
    FH << " limits " << ilower << " " << iupper << endl;
    HYPRE_IJMatrixCreate(HypreComm,ilower,iupper,jlower,jupper,&ij_matrix);
    HYPRE_IJMatrixSetObjectType(ij_matrix,HYPRE_PARCSR);

    int ij_index,udsize,neighs;
    const int total_rows=iupper-ilower+1;
    int* maxcols=new int[total_rows];
    int countr=0;
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p==0)
	  maxcols[countr]=1;
	else
	  {
	    neighs=p->get_ij_neighbors_size();
	    if(neighs == 0)
	      maxcols[countr]=1;
	    else if(neighs == 1)
	      maxcols[countr]=2;
	    else if(neighs == 6)
	      {
		if(p->get_inside())
		  maxcols[countr]=7;
		else
		  maxcols[countr]=1;
	      }	      
	    else
	      assert(0);
	  }
	countr++;
      }
    HYPRE_IJMatrixSetRowSizes(ij_matrix,maxcols);
    delete [] maxcols;
    HYPRE_IJMatrixInitialize(ij_matrix);
    HYPRE_IJVector ij_vector_pot;
    HYPRE_IJVector ij_vector_rho;
    HYPRE_ParVector par_vector_pot;
    HYPRE_ParVector par_vector_rho;
    HYPRE_IJVectorCreate(HypreComm,jlower,jupper,&ij_vector_pot);
    HYPRE_IJVectorCreate(HypreComm,jlower,jupper,&ij_vector_rho);
    HYPRE_IJVectorSetObjectType(ij_vector_pot,HYPRE_PARCSR);
    HYPRE_IJVectorSetObjectType(ij_vector_rho,HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(ij_vector_pot);
    HYPRE_IJVectorInitialize(ij_vector_rho);
    const double pi = 4.0*atan(1.0);
    const int length=frac.get_grid_length();
    double g_c=4.0*pi/static_cast<double>(length*length)*pow(4.0,-level);
    vector <int> ij_ud(6);
    const int nrows=1;
    int rows[1];
    int ncols[1];
    int cols[7];
    double coef1[1]={1.0};
    double coef2[2]={1.0,-1.0};
    double coef7[7]={-6.0,1.0,1.0,1.0,1.0,1.0,1.0};
    double rhov[1];
    double potv[1];
    double rho,pot;
    const int total=mem.ij_counts[HypreRank];
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p==0)
	  {
	    rows[0]=ilower;
	    ncols[0]=1;
	    cols[0]=ilower;
	    potv[0]=1.0;
	    rhov[0]=1.0;
	    udsize=0;
	    HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef1);
	    HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv);
	    HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov);
	    FH << " null point " << endl;
	    continue;
	  }
	else
	  {
	    p->get_hypre_info(ij_index,ij_ud,rho,pot);
	    udsize=ij_ud.size();
	  }
	if(udsize == 0)
	  {
	    rows[0]=ij_index;
	    ncols[0]=1;
	    cols[0]=ij_index;
	    potv[0]=pot;
	    rhov[0]=pot;
	    HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef1);
	    HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv);
	    HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov);
	  }
	else if(udsize == 1)
	  {
	    rows[0]=ij_index;
	    ncols[0]=2;
	    cols[0]=ij_index;
	    cols[1]=ij_ud[0];
	    potv[0]=pot;
	    rhov[0]=0.0;
	    HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef2);
	    HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv);
	    HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov);
	  }
	else if(udsize == 6)
	  {
	    inside=p->get_inside();
	    rows[0]=ij_index;
	    cols[0]=ij_index;
	    potv[0]=pot;
	    if(inside)
	      {
		ncols[0]=7;
		for(int ni=0;ni<6;ni++)
		  cols[ni+1]=ij_ud[ni];
		rhov[0]=rho*g_c;
		HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef7);
	      }
	    else
	      {
		ncols[0]=1;
		rhov[0]=pot;
		HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef1);
	      }
	    HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv);
	    HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov);
	  }
	else
	  assert(0);
      }
    HYPRE_IJMatrixAssemble(ij_matrix);
    HYPRE_IJMatrixGetObject(ij_matrix,(void **) &par_matrix);
    HYPRE_IJVectorAssemble(ij_vector_pot);
    HYPRE_IJVectorAssemble(ij_vector_rho);
    HYPRE_IJVectorGetObject(ij_vector_pot,(void **) &par_vector_pot);
    HYPRE_IJVectorGetObject(ij_vector_rho,(void **) &par_vector_rho);

    HYPRE_Solver par_precond;
    HYPRE_BoomerAMGCreate(&par_precond);
    HYPRE_BoomerAMGSetPrintLevel(par_precond, 1);
    HYPRE_BoomerAMGSetPrintFileName(par_precond, "amg_pre.log");
    HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
    HYPRE_BoomerAMGSetNumSweeps(par_precond, 1);
    HYPRE_BoomerAMGSetRelaxType(par_precond, 6); /* Sym G.S./Jacobi hybrid */ 
    HYPRE_BoomerAMGSetTol(par_precond, 0.0);
    HYPRE_BoomerAMGSetMaxIter(par_precond, 1);


    HYPRE_Solver par_solver;
    HYPRE_ParCSRPCGCreate(HypreComm, &par_solver);

      /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_PCGSetMaxIter(par_solver, frac.get_maxits()); /* max iterations */
    HYPRE_PCGSetTol(par_solver, frac.get_epsilon_sor()); /* conv. tolerance */
    HYPRE_PCGSetTwoNorm(par_solver, 1); /* use the two norm as the stopping criteria */
    HYPRE_PCGSetPrintLevel(par_solver, 2); /* print solve info */
    HYPRE_PCGSetLogging(par_solver, 1); /* needed to get run info later */

    /* Set the PCG preconditioner */
    HYPRE_PCGSetPrecond(par_solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, par_precond);

    /* Now setup and solve! */
    HYPRE_ParCSRPCGSetup(par_solver, par_matrix, par_vector_rho, par_vector_pot);
    HYPRE_ParCSRPCGSolve(par_solver, par_matrix, par_vector_rho, par_vector_pot);

      /* Run info - needed logging turned on */
    int its;
    double final_res_norm;
    HYPRE_PCGGetNumIterations(par_solver, &its);
    HYPRE_PCGGetFinalRelativeResidualNorm(par_solver, &final_res_norm);

    FH << "fini " << level << " " << total << " " << its << " " << final_res_norm << endl;
    assert(its < frac.get_maxits());
    HYPRE_IJMatrixDestroy(ij_matrix);
    HYPRE_IJVectorDestroy(ij_vector_rho);
    HYPRE_BoomerAMGDestroy(par_precond);
    HYPRE_ParCSRPCGDestroy(par_solver);
    double pot0=-1.0;
    int ni=mem.ij_offsets[HypreRank];
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p == 0)
	  FH << " OUT0" << endl;
	else
	  {
	    rows[0]=ni;
	    HYPRE_IJVectorGetValues(ij_vector_pot,1,rows,potv);
	    p->set_potential_point(potv[0]);
	  }
	ni++;
      }
    HYPRE_IJVectorDestroy(ij_vector_pot);
    FH << " exit hypre solver " << level << " " << total << endl;
  }
}
