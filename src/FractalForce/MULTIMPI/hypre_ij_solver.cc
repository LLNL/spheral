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
  void hypre_ij_solver(Fractal& frac,Fractal_Memory& mem,int level,bool& do_over)
  {
    static vector <double> Hypre_sum_time(frac.get_level_max()+1,0.0);
    double Hypre_total_time=0.0;
    double Hypre_search_time=0.0;
    double Hypre_gen_time=0.0;
    double Hypre_setup_time=0.0;
    double Hypre_solve_time=0.0;
    double Hypre_dump_time=0.0;
    const int FractalRank=mem.p_mess->FractalRank;
    if(FractalRank == 0)
      cout << "Hypre Calc " << mem.steps << " " << level << endl;
    ofstream& FH=mem.p_file->FileHypre;
    ofstream& FHT=mem.p_file->FileHypreTime;
    FH << " enter hypre solver " << level << " steps " << mem.steps << endl;
    Hypre_total_time=-mem.p_mess->Clock();
    Hypre_search_time=-mem.p_mess->Clock();
    vector <Point*>hypre_points;
    if(!hypre_ij_numbering(mem,frac,hypre_points,level))
      {
	FH << " nothing here hypre solver " << level << endl;
	return;
      }
    Hypre_search_time+=mem.p_mess->Clock();

    FHT << endl;
    FHT << scientific;
    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_search_time << "\t" << "Search Time" << endl;

    Hypre_gen_time=-mem.p_mess->Clock();
    FH << " really enter hypre solver a " << level << endl;
    bool inside;
    MPI_Comm HypreComm=mem.p_mess->HypreWorld;
    HYPRE_IJMatrix ij_matrix;
    HYPRE_ParCSRMatrix par_matrix;
    const int ilower=mem.ij_offsets[FractalRank];
    const int iupper=ilower+mem.ij_counts[FractalRank]-1;
    const int jlower=ilower;
    const int jupper=iupper;
    FH << " limits " << ilower << " " << iupper << endl;
    hypre_eror(FH,level,0,HYPRE_IJMatrixCreate(HypreComm,ilower,iupper,jlower,jupper,&ij_matrix));
    hypre_eror(FH,level,1,HYPRE_IJMatrixSetObjectType(ij_matrix,HYPRE_PARCSR));

    int ij_index,udsize,neighs;
    const int total_rows=iupper-ilower+1;
    int* maxcols=new int[total_rows];
    int countr=0;
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p)
	  {
	    neighs=p->get_ij_neighbors_size();
	    if(neighs == 0)
	      maxcols[countr]=1;
	    else if(neighs == 1)
	      maxcols[countr]=2;
	    else if(neighs == 6)
	      {
		assert(p->get_inside());
		maxcols[countr]=7;
	      }	      
	    else
	      assert(0);
	  }
	else
	  maxcols[countr]=1;
	countr++;
      }
    assert(countr == total_rows);
    HYPRE_IJMatrixSetRowSizes(ij_matrix,maxcols);
    delete [] maxcols;
    HYPRE_IJMatrixSetMaxOffProcElmts(ij_matrix,0);
    hypre_eror(FH,level,2,HYPRE_IJMatrixInitialize(ij_matrix));
    HYPRE_IJVector ij_vector_pot;
    HYPRE_IJVector ij_vector_rho;
    HYPRE_ParVector par_vector_pot;
    HYPRE_ParVector par_vector_rho;
    FH << " really enter hypre solver b " << level << endl;
    hypre_eror(FH,level,3,HYPRE_IJVectorCreate(HypreComm,jlower,jupper,&ij_vector_pot));
    hypre_eror(FH,level,4,HYPRE_IJVectorCreate(HypreComm,jlower,jupper,&ij_vector_rho));
    hypre_eror(FH,level,5,HYPRE_IJVectorSetObjectType(ij_vector_pot,HYPRE_PARCSR));
    hypre_eror(FH,level,6,HYPRE_IJVectorSetObjectType(ij_vector_rho,HYPRE_PARCSR));
    hypre_eror(FH,level,7,HYPRE_IJVectorInitialize(ij_vector_pot));
    hypre_eror(FH,level,8,HYPRE_IJVectorInitialize(ij_vector_rho));
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
	    hypre_eror(FH,level,9,HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef1));
	    hypre_eror(FH,level,10,HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv));
	    hypre_eror(FH,level,11,HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov));
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
	    hypre_eror(FH,level,12,HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef1));
	    hypre_eror(FH,level,13,HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv));
	    hypre_eror(FH,level,14,HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov));
	  }
	else if(udsize == 1)
	  {
	    rows[0]=ij_index;
	    ncols[0]=2;
	    cols[0]=ij_index;
	    cols[1]=ij_ud[0];
	    potv[0]=pot;
	    rhov[0]=0.0;
	    hypre_eror(FH,level,15,HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef2));
	    hypre_eror(FH,level,16,HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv));
	    hypre_eror(FH,level,17,HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov));
	  }
	else if(udsize == 6)
	  {
	    inside=p->get_inside();
	    rows[0]=ij_index;
	    cols[0]=ij_index;
	    potv[0]=pot;
	    assert(inside);
	    ncols[0]=7;
	    for(int ni=0;ni<6;ni++)
	      cols[ni+1]=ij_ud[ni];
	    rhov[0]=rho*g_c;
	    hypre_eror(FH,level,18,HYPRE_IJMatrixSetValues(ij_matrix,nrows,ncols,rows,cols,coef7));
	    hypre_eror(FH,level,20,HYPRE_IJVectorSetValues(ij_vector_pot,1,rows,potv));
	    hypre_eror(FH,level,21,HYPRE_IJVectorSetValues(ij_vector_rho,1,rows,rhov));
	  }
	else
	  assert(0);
      }
    hypre_eror(FH,level,22,HYPRE_IJMatrixAssemble(ij_matrix));
    hypre_eror(FH,level,23,HYPRE_IJMatrixGetObject(ij_matrix,(void **) &par_matrix));
    hypre_eror(FH,level,24,HYPRE_IJVectorAssemble(ij_vector_pot));
    hypre_eror(FH,level,25,HYPRE_IJVectorAssemble(ij_vector_rho));
    hypre_eror(FH,level,26,HYPRE_IJVectorGetObject(ij_vector_pot,(void **) &par_vector_pot));
    hypre_eror(FH,level,27,HYPRE_IJVectorGetObject(ij_vector_rho,(void **) &par_vector_rho));
    HYPRE_Solver par_solver;
    hypre_eror(FH,level,28,HYPRE_BoomerAMGCreate(&par_solver));
    hypre_eror(FH,level,-1,HYPRE_BoomerAMGSetDebugFlag(par_solver,1));
    hypre_eror(FH,level,29,HYPRE_BoomerAMGSetCoarsenType(par_solver, 6));
    hypre_eror(FH,level,30,HYPRE_BoomerAMGSetStrongThreshold(par_solver, 0.55));
    hypre_eror(FH,level,31,HYPRE_BoomerAMGSetTol(par_solver, frac.get_epsilon_sor()));
    hypre_eror(FH,level,32,HYPRE_BoomerAMGSetPrintLevel(par_solver, 1));
    hypre_eror(FH,level,33,HYPRE_BoomerAMGSetPrintFileName(par_solver, "amg_real.log"));
    hypre_eror(FH,level,34,HYPRE_BoomerAMGSetMaxIter(par_solver, frac.get_maxits()));
    Hypre_gen_time+=mem.p_mess->Clock();

    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_gen_time << "\t" << "Gen    Time" << endl;

    Hypre_setup_time=-mem.p_mess->Clock();
    hypre_eror(FH,level,35,HYPRE_BoomerAMGSetup(par_solver, par_matrix, par_vector_rho, par_vector_pot));
    Hypre_setup_time+=mem.p_mess->Clock();

    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_setup_time << "\t" << "Setup  Time" << endl;

    Hypre_solve_time=-mem.p_mess->Clock();
    hypre_eror(FH,level,36,HYPRE_BoomerAMGSolve(par_solver, par_matrix, par_vector_rho, par_vector_pot));
    Hypre_solve_time+=mem.p_mess->Clock();

    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_solve_time << "\t" << "Solve  Time" << endl;

    Hypre_dump_time=-mem.p_mess->Clock();
    int its;
    double final_res_norm;
    hypre_eror(FH,level,37,HYPRE_BoomerAMGGetNumIterations(par_solver, &its));
    hypre_eror(FH,level,38,HYPRE_BoomerAMGGetFinalRelativeResidualNorm(par_solver,&final_res_norm));
    FH << "fini " << level << " " << total << " " << its << " " << final_res_norm << endl;

    do_over=its >= frac.get_maxits();
    if(do_over)
      FHT << " no convergence, try again " << " " << level << endl;
    hypre_eror(FH,level,39,HYPRE_IJMatrixDestroy(ij_matrix));
    hypre_eror(FH,level,40,HYPRE_IJVectorDestroy(ij_vector_rho));
    hypre_eror(FH,level,41,HYPRE_BoomerAMGDestroy(par_solver));

    if(!do_over)
      {
	int ni=mem.ij_offsets[FractalRank];
	for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
	  {
	    Point* p=*point_itr;
	    if(p)
	      {
		rows[0]=ni;
		hypre_eror(FH,level,42,HYPRE_IJVectorGetValues(ij_vector_pot,1,rows,potv));
		p->set_potential_point(potv[0]);
	      }
	    else
	      FH << " OUT0" << endl;
	    ni++;
	  }
      }
    hypre_eror(FH,level,43,HYPRE_IJVectorDestroy(ij_vector_pot));
    Hypre_dump_time+=mem.p_mess->Clock();

    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_dump_time << "\t" << "Dump   Time" << endl;

    Hypre_total_time+=mem.p_mess->Clock();
    Hypre_sum_time[level]+=Hypre_total_time;

    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_total_time << "\t" << Hypre_sum_time[level] << " Total Time " << endl;

    FH << " exit hypre solver " << level << " " << total << " steps " << mem.steps << endl;
  }
  void hypre_eror(ofstream& FH,int level,int ni,int er)
  {
    if(er == 0)
      return;
    char describe[100];
    FH << " hypre describe " << ni << " " << level << " " << er << " ";
    HYPRE_DescribeError(er,describe);
    FH << describe << endl;
  }
}
