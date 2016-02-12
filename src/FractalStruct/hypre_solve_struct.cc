#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#ifndef _Hypre_Defined_
#define _Hypre_Defined_
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE_struct_ls.h"
#endif
namespace FractalSpace
{
  void hypre_solve_struct(Fractal_Memory& mem,int level,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   Amatrix;
    HYPRE_StructVector   rho;
    HYPRE_StructVector   pot;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;
    vector <int>BBox=mem.FRBBoxesLev[level];
    double pi = 4.0*atan(1.0);
    int length=mem.p_fractal->get_grid_length();
    double g_c=4.0*pi/static_cast<double>(length*length)*pow(4.0,-level);
    int spacing=Misc::pow(2,mem.p_fractal->get_level_max()-level);
    int lowerBOX[3];
    int upperBOX[3];
    HYPRE_StructGridCreate(mem.p_mess->HypreWorld,3,&grid);
    for(int B=0;B<SBoxes.size();B++)
      {
	int ni=0;
	for(int ni2=0;ni2<6;ni2+=2)
	  {
	    lowerBOX[ni]=SBoxes[B][ni2]/spacing;
	    upperBOX[ni]=SBoxes[B][ni2+1]/spacing-1;
	    ni++;
	  }
	HYPRE_StructGridSetExtents(grid,lowerBOX,upperBOX);
      }
    HYPRE_StructGridAssemble(grid);
    HYPRE_StructStencilCreate(3,7,&stencil);
    int offsets[7][3] = {{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    for (int entry = 0; entry < 7; entry++)
      HYPRE_StructStencilSetElement(stencil,entry,offsets[entry]);
    HYPRE_StructMatrixCreate(mem.p_mess->HypreWorld,grid,stencil, &Amatrix);
    HYPRE_StructMatrixInitialize(Amatrix);
    int stencil_indices[7] = {0,1,2,3,4,5,6};
    for(int B=0;B<SBoxes.size();B++)
      {
	vector <int>pos(3);
	vector <double>values;
	int VOL=1;
	int ni=0;
	for(int ni2=0;ni2<6;ni2+=2)
	  {
	    lowerBOX[ni]=(SBoxes[B][ni2]/spacing);
	    upperBOX[ni]=(SBoxes[B][ni2+1]/spacing-1);
	    VOL*=upperBOX[ni]-lowerBOX[ni]+1;
	    ni++;
	  }
	vector <Point*>::iterator p_itr=SPoints[B].begin();
	for (int i = 0; i < VOL;i++)
	  {
	    Point* p=*p_itr;
	    values.push_back(-6.0);
	    for (int j = 0; j < 6; j++)
	      {
		Point* p1=p->get_point_ud(j);
		p1->get_pos_point(pos);
		if(p1->get_inside() || pos[j/2] == BBox[j])
		  values.push_back(1.0);
		else
		  values.push_back(0.0);
	      }
	    p_itr++;
	  }
	HYPRE_StructMatrixAddToBoxValues(Amatrix,lowerBOX,upperBOX,7,
				       stencil_indices,&(*values.begin()));
	values.clear();
      }
    HYPRE_StructMatrixAssemble(Amatrix);
    HYPRE_StructVectorCreate(mem.p_mess->HypreWorld, grid, &rho);
    HYPRE_StructVectorCreate(mem.p_mess->HypreWorld, grid, &pot);
    
    HYPRE_StructVectorInitialize(rho);
    HYPRE_StructVectorInitialize(pot);
    for(int B=0;B<SBoxes.size();B++)
      {
	vector <double>dens_values;
	vector <double>pot_values;
	int VOL=1;
	int ni=0;
	for(int ni2=0;ni2<6;ni2+=2)
	  {
	    lowerBOX[ni]=SBoxes[B][ni2]/spacing;
	    upperBOX[ni]=SBoxes[B][ni2+1]/spacing-1;
	    VOL*=upperBOX[ni]-lowerBOX[ni]+1;
	    ni++;
	  }
	vector <Point*>::iterator p_itr=SPoints[B].begin();
	for (int i = 0; i < VOL; i++)
	  {
	    vector <int>pos(3);
	    Point* p=*p_itr;
	    double density=p->get_density_point()*g_c;
	    for(int ni=0;ni<6;ni++)
	      {
		Point* p1=p->get_point_ud(ni);
		if(p1->get_inside() || pos[ni/2] == BBox[ni])
		  continue;
		density-=p1->get_potential_point();
	      }
	    dens_values.push_back(density);
	    pot_values.push_back(p->get_potential_point());
	    p_itr++;
	  }
	HYPRE_StructVectorAddToBoxValues(rho,lowerBOX,upperBOX,&(*dens_values.begin()));
	HYPRE_StructVectorAddToBoxValues(pot,lowerBOX,upperBOX,&(*pot_values.begin()));
	dens_values.clear();
	pot_values.clear();
      }
    HYPRE_StructVectorAssemble(rho);
    HYPRE_StructVectorAssemble(pot);
    
    HYPRE_StructPCGCreate(mem.p_mess->HypreWorld, &solver);
    HYPRE_StructPCGSetMaxIter(solver, 200 );
    HYPRE_StructPCGSetTol(solver, 1.0e-06 );
    HYPRE_StructPCGSetTwoNorm(solver, 1 );
    HYPRE_StructPCGSetRelChange(solver, 0 );
    HYPRE_StructPCGSetPrintLevel(solver, 2 );
  
    HYPRE_StructPFMGCreate(mem.p_mess->HypreWorld, &precond);
    HYPRE_StructPFMGSetMaxIter(precond, 1);
    HYPRE_StructPFMGSetTol(precond, 0.0);
    //    HYPRE_StructPFMGSetZeroGuess(precond);
    HYPRE_StructPFMGSetRAPType(precond, 1);
    //    HYPRE_StructPFMGSetRelaxType(precond, relax);
    //    HYPRE_StructPFMGSetNumPreRelax(precond, n_pre);
    //    HYPRE_StructPFMGSetNumPostRelax(precond, n_post);
    HYPRE_StructPFMGSetSkipRelax(precond, 1);
    HYPRE_StructPFMGSetPrintLevel(precond, 0);
    HYPRE_StructPFMGSetLogging(precond, 0);
    HYPRE_StructPCGSetPrecond(solver,
			      HYPRE_StructPFMGSolve,
			      HYPRE_StructPFMGSetup,
			      precond);
    HYPRE_StructPCGSetup(solver,Amatrix,rho,pot);
    HYPRE_StructPCGSolve(solver,Amatrix,rho,pot);
    int num_iterations=-1;
    double final_res_norm=-1.0;
    HYPRE_StructPCGGetNumIterations(solver,&num_iterations );
    HYPRE_StructPCGGetFinalRelativeResidualNorm( solver, &final_res_norm );
    HYPRE_StructPCGDestroy(solver);
    HYPRE_StructPFMGDestroy(precond);
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructMatrixDestroy(Amatrix);
    HYPRE_StructVectorDestroy(rho);
    for(int B=0;B<SBoxes.size();B++)
      {
	vector <double>pot_values;
	int ni=0;
	int VOL=1;
	for(int ni2=0;ni2<6;ni2+=2)
	  {
	    lowerBOX[ni]=SBoxes[B][ni2]/spacing;
	    upperBOX[ni]=SBoxes[B][ni2+1]/spacing-1;
	    VOL*=upperBOX[ni]-lowerBOX[ni]+1;
	    ni++;
	  }
	pot_values.resize(VOL);
	HYPRE_StructVectorGetBoxValues(pot,lowerBOX,upperBOX,&(*pot_values.begin()));
	vector <Point*>::iterator p_itr=SPoints[B].begin();
	for (int i = 0; i < VOL; i++)
	  {
	    (*p_itr)->set_potential_point(pot_values[i]);
	    p_itr++;
	  }
	pot_values.clear();
      }
    //add_buffer_values
  //
    HYPRE_StructVectorDestroy(pot);
  }
}
