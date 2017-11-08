#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
#ifndef _Hypre_Defined_
#define _Hypre_Defined_
#include "HYPRE.h"
#include "_hypre_utilities.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#endif
namespace FractalSpace
{
  void hypre_solver(Fractal& frac,Fractal_Memory& mem,const int& level)
  {
    ofstream& FH=mem.p_file->FileHypre;
    char describe[100];
    vector <int>hyp_error(50,0);
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
    MPI_Comm HypreComm=mem.p_mess->HypreWorld;;
    HYPRE_IJMatrix ij_matrix;
    HYPRE_ParCSRMatrix par_matrix;
    const int ilower=mem.ij_offsets[FractalRank];
    const int iupper=ilower+mem.ij_counts[FractalRank]-1;
    const int jlower=ilower;
    const int jupper=iupper;
    FH << " limits " << ilower << " " << iupper << endl;
    hyp_error[0]=HYPRE_IJMatrixCreate(HypreComm,ilower,iupper,jlower,jupper,&ij_matrix);
    hyp_error[1]=HYPRE_IJMatrixSetObjectType(ij_matrix,HYPRE_PARCSR);

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
    hyp_error[2]=HYPRE_IJMatrixInitialize(ij_matrix);
    HYPRE_IJVector ij_vector_pot;
    HYPRE_IJVector ij_vector_rho;
    HYPRE_ParVector par_vector_pot;
    HYPRE_ParVector par_vector_rho;
    FH << " really enter hypre solver b " << level << endl;
    hyp_error[3]=HYPRE_IJVectorCreate(HypreComm,jlower,jupper,&ij_vector_pot);
    hyp_error[4]=HYPRE_IJVectorCreate(HypreComm,jlower,jupper,&ij_vector_rho);
    hyp_error[5]=HYPRE_IJVectorSetObjectType(ij_vector_pot,HYPRE_PARCSR);
    hyp_error[6]=HYPRE_IJVectorSetObjectType(ij_vector_rho,HYPRE_PARCSR);
    hyp_error[7]=HYPRE_IJVectorInitialize(ij_vector_pot);
    hyp_error[8]=HYPRE_IJVectorInitialize(ij_vector_rho);
    FH << " really enter hypre solver c " << level << endl;
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
    const int total=mem.ij_counts[FractalRank];
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
	    hyp_error[9]=HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef1);
	    hyp_error[10]=HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv);
	    hyp_error[11]=HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov);
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
	    hyp_error[9]=HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef1);
	    hyp_error[10]=HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv);
	    hyp_error[11]=HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov);
	    //	    FH << " UD0 ";
	  }
	else if(udsize == 1)
	  {
	    rows[0]=ij_index;
	    ncols[0]=2;
	    cols[0]=ij_index;
	    cols[1]=ij_ud[0];
	    potv[0]=pot;
	    //	    potv[0]=0.0;
	    rhov[0]=0.0;
	    hyp_error[12]=HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef2);
	    hyp_error[13]=HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv);
	    hyp_error[14]=HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov);
	    //	    FH << " UD1 ";
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
		hyp_error[15]=HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef7);
	      }
	    else
	      {
		ncols[0]=1;
		rhov[0]=pot;
		hyp_error[15]=HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef1);
	      }
	    hyp_error[16]=HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv);
	    hyp_error[17]=HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov);
	    //	    FH << " UD6 ";
	  }
	else
	  assert(0);
	//	for(int ni=0;ni<=udsize;ni++)
	//	  FH << cols[ni] << " ";
	//	FH << potv[0] << " " << rhov[0] << endl;
      }
    FH << " made it this far a " << endl;
    hyp_error[18]=HYPRE_IJMatrixAssemble(ij_matrix);
    FH << " made it this far b " << endl;
    //    HYPRE_IJMatrixPrint(ij_matrix,"matrix.a");
    hyp_error[19]=HYPRE_IJMatrixGetObject(ij_matrix,(void **) &par_matrix);
    FH << " made it this far c " << endl;
    hyp_error[20]=HYPRE_IJVectorAssemble(ij_vector_pot);
    FH << " made it this far d " << endl;
    //    HYPRE_IJVectorPrint(ij_vector_pot,"vector.pot");
    hyp_error[21]=HYPRE_IJVectorAssemble(ij_vector_rho);
    FH << " made it this far e " << endl;
    //    HYPRE_IJVectorPrint(ij_vector_rho,"vector.rho");
    hyp_error[22]=HYPRE_IJVectorGetObject(ij_vector_pot,(void **) &par_vector_pot);
    FH << " made it this far f " << endl;
    hyp_error[23]=HYPRE_IJVectorGetObject(ij_vector_rho,(void **) &par_vector_rho);
    FH << " made it this far g " << endl;
    /*
    HYPRE_Solver par_precond;
    hyp_error[24]=HYPRE_BoomerAMGCreate(&par_precond);
    hyp_error[25]=HYPRE_BoomerAMGSetCoarsenType(par_precond, 6);
    hyp_error[26]=HYPRE_BoomerAMGSetStrongThreshold(par_precond, 0.25);
    hyp_error[27]=HYPRE_BoomerAMGSetTol(par_precond, 0.0);
    hyp_error[28]=HYPRE_BoomerAMGSetPrintLevel(par_precond, 1);
    hyp_error[29]=HYPRE_BoomerAMGSetPrintFileName(par_precond, "amg_pre.log");
    hyp_error[30]=HYPRE_BoomerAMGSetMaxIter(par_precond, 1);
    */
    HYPRE_Solver par_solver;
    hyp_error[31]=HYPRE_BoomerAMGCreate(&par_solver);
    FH << " made it this far h " << endl;
    hyp_error[32]=HYPRE_BoomerAMGSetCoarsenType(par_solver, 6);
    FH << " made it this far i " << endl;
    hyp_error[33]=HYPRE_BoomerAMGSetStrongThreshold(par_solver, 0.55);
    FH << " made it this far j " << endl;
    hyp_error[34]=HYPRE_BoomerAMGSetTol(par_solver, frac.get_epsilon_sor());
    FH << " made it this far k " << endl;
    hyp_error[35]=HYPRE_BoomerAMGSetPrintLevel(par_solver, 1);
    FH << " made it this far l " << endl;
    hyp_error[36]=HYPRE_BoomerAMGSetPrintFileName(par_solver, "amg_real.log");
    FH << " made it this far m " << endl;
    hyp_error[37]=HYPRE_BoomerAMGSetMaxIter(par_solver, frac.get_maxits());
    FH << " made it this far n " << endl;
    /*
    hyp_error[38]=HYPRE_BoomerAMGSetup(par_precond, par_matrix, par_vector_rho, par_vector_pot);
    */
    hyp_error[39]=HYPRE_BoomerAMGSetup(par_solver, par_matrix, par_vector_rho, par_vector_pot);
    FH << " made it this far o " << endl;
    /*
    hyp_error[40]=HYPRE_BoomerAMGSolve(par_precond, par_matrix, par_vector_rho, par_vector_pot);
    */
    hyp_error[41]=HYPRE_BoomerAMGSolve(par_solver, par_matrix, par_vector_rho, par_vector_pot);
    FH << " made it this far p " << endl;

    int its;
    double final_res_norm;
    hyp_error[42]=HYPRE_BoomerAMGGetNumIterations(par_solver, &its);
    FH << " made it this far q " << endl;
    hyp_error[43]=HYPRE_BoomerAMGGetFinalRelativeResidualNorm(par_solver,&final_res_norm);
    FH << "fini " << level << " " << total << " " << its << " " << final_res_norm << endl;
    assert(its < frac.get_maxits());
    hyp_error[44]=HYPRE_IJMatrixDestroy(ij_matrix);
    hyp_error[45]=HYPRE_IJVectorDestroy(ij_vector_rho);
    /*
    hyp_error[46]=HYPRE_BoomerAMGDestroy(par_precond);
    */
    hyp_error[47]=HYPRE_BoomerAMGDestroy(par_solver);
    double pot0=-1.0;
    int ni=mem.ij_offsets[FractalRank];
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p == 0)
	  FH << " OUT0" << endl;
	else
	  {
	    rows[0]=ni;
	    hyp_error[48]=HYPRE_IJVectorGetValues(ij_vector_pot,1,rows,potv);
	    //	    pot0=p->get_potential_point();
	    //	    inside=p->get_inside();
	    //	    edge=p->get_edge_point();
	    //	    buff=p->get_buffer_point();
	    //	    pass=p->get_passive_point();
	    p->set_potential_point(potv[0]);
	    //	    FH << " OUT " << ni << " " << inside << edge << buff << pass << " " << pot0 << " " << potv[0]-pot0 << endl;
	  }
	ni++;
      }
    hyp_error[49]=HYPRE_IJVectorDestroy(ij_vector_pot);
    for(int ni=0;ni<50;ni++)
      if(hyp_error[ni]) hypre_eror(FH,hyp_error[ni],level,ni);
    FH << " exit hypre solver " << level << " " << total << endl;
  }
  void hypre_eror(ofstream& FH,const int& er,const int& level,const int& lab)
  {
    char describe[100];
    FH << " hypre error " << lab << " " << level << " " << er << endl;
    HYPRE_DescribeError(er,describe);
    FH << describe << endl;
  }
}
