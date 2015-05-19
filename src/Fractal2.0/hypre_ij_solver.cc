#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_ij_solver(Fractal& frac,Fractal_Memory& mem,const int& level)
  {
    static vector <double> Hypre_sum_time(frac.get_level_max()+1,0.0);
    double Hypre_total_time=0.0;
    double Hypre_search_time=0.0;
    double Hypre_gen_time=0.0;
    double Hypre_setup_time=0.0;
    double Hypre_solve_time=0.0;
    double Hypre_dump_time=0.0;
    mem.p_mess->HypreRank=-1;
    mem.p_mess->IAmAHypreNode=false;;
    int FractalRank=mem.p_mess->FractalRank;
    int HypreRank=-1;
    if(FractalRank ==0)
      cerr << "Hypre Calc " << mem.steps << " " << level << " " << FractalRank << " " << sizeof(HYPRE_Int) << " Hypre_Int" << "\n";
    FILE* PFH=mem.p_file->PFHypre;
    ofstream& FHT=mem.p_file->DUMPS;
    fprintf(PFH," enter hypre solver %d %d \n",level,mem.steps);
    Hypre_total_time=-mem.p_mess->Clock();
    Hypre_search_time=-mem.p_mess->Clock();
    vector <Point*>hypre_points;
    frac.timing(-1,32);
    if(!hypre_ij_numbering(mem,hypre_points,level))
      {
	//	fprintf(PFH," nothing here hypre solver %d %d %d \n",level,mem.p_mess->HypreRank,FractalRank);
	frac.timing(1,32);
	return;
      }
    frac.timing(1,32);
    FHT << " TAKING STEPS " << mem.steps << " " << level << " " << FractalRank << " " << mem.p_mess->HypreRank;
    if(!mem.Touchy.empty())
      FHT << " " << mem.Touchy.front();
    FHT << "\n";
    Hypre_search_time+=mem.p_mess->Clock();
    HypreRank=mem.p_mess->HypreRank;
    if(mem.p_mess->IAmAHypreNode)
      Full_Stop(mem,mem.p_mess->HypreWorld,37);
    FHT << "\n";
    FHT << scientific;
    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_search_time << "\t" << "Search Time Buffer" << "\n";
    int total=0;
    //    fprintf(PFH," Am I a Hypre Node %d %d %d \n",mem.p_mess->IAmAHypreNode,HypreRank,FractalRank);
    if(mem.p_mess->IAmAHypreNode)
      {
	Hypre_gen_time=-mem.p_mess->Clock();
	//	fprintf(PFH," really enter hypre solver a %d \n",level);
	bool load_balance=false;
	HYPRE_Int off_elements=hypre_load_balance(mem,hypre_points,load_balance);
	MPI_Comm HypreComm=mem.p_mess->HypreWorld;
	HYPRE_IJMatrix ij_matrix;
	HYPRE_ParCSRMatrix par_matrix;
	const HYPRE_Int ilower=mem.ij_offsetsB[HypreRank];
	const HYPRE_Int iupper=mem.ij_offsetsB[HypreRank+1]-1;
	const HYPRE_Int jlower=ilower;
	const HYPRE_Int jupper=iupper;
	fprintf(PFH," limits %d %d \n",ilower,iupper);
	hypre_eror(PFH,level,0,HYPRE_IJMatrixCreate(HypreComm,ilower,iupper,jlower,jupper,&ij_matrix));
	hypre_eror(PFH,level,1,HYPRE_IJMatrixSetObjectType(ij_matrix,HYPRE_PARCSR));
	HYPRE_IJMatrixSetMaxOffProcElmts(ij_matrix,off_elements);
	HYPRE_Int ij_index,udsize;
	vector <HYPRE_Int> maxcols(mem.ij_countsB[HypreRank],7);
	hypre_eror(PFH,level,-1,HYPRE_IJMatrixSetRowSizes(ij_matrix,&(*maxcols.begin())));
	maxcols.clear();
	hypre_eror(PFH,level,2,HYPRE_IJMatrixInitialize(ij_matrix));
	//	fprintf(PFH," really enter hypre solver a %d \n",level);
	HYPRE_IJVector ij_vector_pot;
	HYPRE_IJVector ij_vector_rho;
	HYPRE_ParVector par_vector_pot;
	HYPRE_ParVector par_vector_rho;
	//	fprintf(PFH," really enter hypre solver b %d \n",level);
	hypre_eror(PFH,level,3,HYPRE_IJVectorCreate(HypreComm,jlower,jupper,&ij_vector_pot));
	hypre_eror(PFH,level,4,HYPRE_IJVectorCreate(HypreComm,jlower,jupper,&ij_vector_rho));
	hypre_eror(PFH,level,5,HYPRE_IJVectorSetObjectType(ij_vector_pot,HYPRE_PARCSR));
	hypre_eror(PFH,level,6,HYPRE_IJVectorSetObjectType(ij_vector_rho,HYPRE_PARCSR));
	hypre_eror(PFH,level,-6,HYPRE_IJVectorSetMaxOffProcElmts(ij_vector_pot,off_elements));
	hypre_eror(PFH,level,-7,HYPRE_IJVectorSetMaxOffProcElmts(ij_vector_rho,off_elements));
	hypre_eror(PFH,level,7,HYPRE_IJVectorInitialize(ij_vector_pot));
	hypre_eror(PFH,level,8,HYPRE_IJVectorInitialize(ij_vector_rho));
	//	fprintf(PFH," really enter hypre solver c %d \n",level);
	const double pi = 4.0*atan(1.0);
	const int length=frac.get_grid_length();
	double g_c=4.0*pi/static_cast<double>(length*length)*pow(4.0,-level);
	vector <HYPRE_Int> ij_ud(6);
	const HYPRE_Int nrows=1;
	HYPRE_Int rows[1];
	HYPRE_Int ncols[1];
	HYPRE_Int cols[7];
	double coef1[1]={1.0};
	double coef2[2]={1.0,-1.0};
	double coef7[7]={-6.0,1.0,1.0,1.0,1.0,1.0,1.0};
	double rhov[1];
	double potv[1];
	double rho,pot;
	total=hypre_points.size();
	for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
	  {
	    vector <Point*> ud6(6);
	    Point* p=*point_itr;
	    p->get_hypre_info(ij_index,ij_ud,rho,pot);
	    udsize=ij_ud.size();
	    if(udsize == 0)
	      {
		rows[0]=ij_index;
		ncols[0]=1;
		cols[0]=ij_index;
		potv[0]=pot;
		rhov[0]=pot;
		hypre_eror(PFH,level,12,HYPRE_IJMatrixAddToValues(ij_matrix,nrows,ncols,rows,cols,coef1));
		hypre_eror(PFH,level,13,HYPRE_IJVectorAddToValues(ij_vector_pot,1,rows,potv));
		hypre_eror(PFH,level,14,HYPRE_IJVectorAddToValues(ij_vector_rho,1,rows,rhov));
	      }
	    else if(udsize == 1)
	      {
		rows[0]=ij_index;
		ncols[0]=2;
		cols[0]=ij_index;
		cols[1]=ij_ud[0];
		potv[0]=pot;
		rhov[0]=0.0;
		hypre_eror(PFH,level,15,HYPRE_IJMatrixAddToValues(ij_matrix,nrows,ncols,rows,cols,coef2));
		hypre_eror(PFH,level,16,HYPRE_IJVectorAddToValues(ij_vector_pot,1,rows,potv));
		hypre_eror(PFH,level,17,HYPRE_IJVectorAddToValues(ij_vector_rho,1,rows,rhov));
	      }
	    else if(udsize == 6)
	      {
		p->get_point_ud(ud6);
		assert(p->get_inside());
		rows[0]=ij_index;
		cols[0]=ij_index;
		potv[0]=pot;
		rhov[0]=rho*g_c;
		HYPRE_Int neighs=1;
		for(int ni=0;ni<6;ni++)
		  {
		    if(ud6[ni]->get_inside())
		      {
			cols[neighs]=ij_ud[ni];
			neighs++;
		      }
		    else
		      rhov[0]-=ud6[ni]->get_potential_point();
		  }
		ncols[0]=neighs;
		hypre_eror(PFH,level,18,HYPRE_IJMatrixAddToValues(ij_matrix,nrows,ncols,rows,cols,coef7));
		hypre_eror(PFH,level,20,HYPRE_IJVectorAddToValues(ij_vector_pot,1,rows,potv));
		hypre_eror(PFH,level,21,HYPRE_IJVectorAddToValues(ij_vector_rho,1,rows,rhov));
	      }
	    else 
	      assert(0);
	  }
	hypre_eror(PFH,level,22,HYPRE_IJMatrixAssemble(ij_matrix));
	hypre_eror(PFH,level,23,HYPRE_IJMatrixGetObject(ij_matrix,(void **) &par_matrix));
	hypre_eror(PFH,level,24,HYPRE_IJVectorAssemble(ij_vector_pot));
	hypre_eror(PFH,level,25,HYPRE_IJVectorAssemble(ij_vector_rho));
	hypre_eror(PFH,level,26,HYPRE_IJVectorGetObject(ij_vector_pot,(void **) &par_vector_pot));
	hypre_eror(PFH,level,27,HYPRE_IJVectorGetObject(ij_vector_rho,(void **) &par_vector_rho));
	HYPRE_Solver par_solver;
	hypre_eror(PFH,level,28,HYPRE_BoomerAMGCreate(&par_solver));
	hypre_eror(PFH,level,-1,HYPRE_BoomerAMGSetDebugFlag(par_solver,0));
	hypre_eror(PFH,level,29,HYPRE_BoomerAMGSetCoarsenType(par_solver, 6));
	hypre_eror(PFH,level,30,HYPRE_BoomerAMGSetStrongThreshold(par_solver, 0.55));
	hypre_eror(PFH,level,31,HYPRE_BoomerAMGSetTol(par_solver, frac.get_epsilon_sor()));
	hypre_eror(PFH,level,32,HYPRE_BoomerAMGSetPrintLevel(par_solver, 0));
	hypre_eror(PFH,level,33,HYPRE_BoomerAMGSetPrintFileName(par_solver, "amg_real.log"));
	hypre_eror(PFH,level,34,HYPRE_BoomerAMGSetMaxIter(par_solver, frac.get_maxits()));
	Hypre_gen_time+=mem.p_mess->Clock();
	FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_gen_time << "\t" << "Gen    Time Buffer" << "\n";

	Hypre_setup_time=-mem.p_mess->Clock();
	hypre_eror(PFH,level,35,HYPRE_BoomerAMGSetup(par_solver, par_matrix, par_vector_rho, par_vector_pot));
	Hypre_setup_time+=mem.p_mess->Clock();

	FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_setup_time << "\t" << "Setup  Time Buffer" << "\n";

	Hypre_solve_time=-mem.p_mess->Clock();
	hypre_eror(PFH,level,36,HYPRE_BoomerAMGSolve(par_solver, par_matrix, par_vector_rho, par_vector_pot));
	Hypre_solve_time+=mem.p_mess->Clock();

	FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_solve_time << "\t" << "Solve  Time Buffer" << "\n";

	Hypre_dump_time=-mem.p_mess->Clock();
	HYPRE_Int its;
	double final_res_norm;
	hypre_eror(PFH,level,37,HYPRE_BoomerAMGGetNumIterations(par_solver, &its));
	hypre_eror(PFH,level,38,HYPRE_BoomerAMGGetFinalRelativeResidualNorm(par_solver,&final_res_norm));
	fprintf(PFH,"fini %d %d %d %12.4E \n",level,total,its,final_res_norm);
	hypre_eror(PFH,level,39,HYPRE_IJMatrixDestroy(ij_matrix));
	hypre_eror(PFH,level,40,HYPRE_IJVectorDestroy(ij_vector_rho));
	hypre_eror(PFH,level,41,HYPRE_BoomerAMGDestroy(par_solver));

	bool do_over=its >= frac.get_maxits();
	if(do_over)
	  {
	    HYPRE_IJMatrixPrint(ij_matrix,"ij_matrix_dump");
	    HYPRE_IJVectorPrint(ij_vector_rho,"ij_vector_rho_dump");
	    HYPRE_IJVectorPrint(ij_vector_pot,"ij_vector_pot_dump");
	    FHT << " no convergence, try again " << " " << level << "\n";
	    assert(its < frac.get_maxits());
	  }
	vector <double>potH(mem.ij_countsB[HypreRank]);
	vector <HYPRE_Int>rowsH(mem.ij_countsB[HypreRank]);
	HYPRE_Int holy_grail=mem.ij_offsetsB[HypreRank];
	for(HYPRE_Int ni=0;ni<mem.ij_countsB[HypreRank];ni++)
	  rowsH[ni]=ni+holy_grail;
	hypre_eror(PFH,level,42,HYPRE_IJVectorGetValues(ij_vector_pot,mem.ij_countsB[HypreRank],&(*rowsH.begin()),&(*potH.begin())));
	rowsH.clear();
	if(load_balance)
	  hypre_send_pots(mem,hypre_points,potH);
	else
	  {
	    HYPRE_Int Mr_smokes_too_much=0;
	    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
	      {
		Point* p=*point_itr;
		rows[0]=p->get_ij_number();
		hypre_eror(PFH,level,42,HYPRE_IJVectorGetValues(ij_vector_pot,1,rows,potv));
		p->set_potential_point(potH[Mr_smokes_too_much]);
		Mr_smokes_too_much++;
	      }
	  }
	potH.clear();
	hypre_eror(PFH,level,43,HYPRE_IJVectorDestroy(ij_vector_pot));
	Hypre_dump_time+=mem.p_mess->Clock();

	FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_dump_time << "\t" << "Dump   Time Buffer" << "\n";
      }
    Hypre_total_time+=mem.p_mess->Clock();
    Hypre_sum_time[level]+=Hypre_total_time;

    FHT << " S" << mem.steps << "S " << "L" << level << "L" << "\t" << Hypre_total_time << "\t" << Hypre_sum_time[level] << " Total Time Buffer ";
    FHT << "\t" <<mem.ij_offsets[mem.p_mess->HypreNodes] << "\n";

    //    fprintf(PFH," exit hypre solver %d %d steps %d \n",level,total, mem.steps);
    mem.p_mess->HypreGroupFree();
  }
  void hypre_eror(FILE* PFH,int level,int ni,int er)
  {
    if(er == 0)
      return;
    char describe[100];
    fprintf(PFH," hypre describe %d %d %d",ni,level,er);
    HYPRE_DescribeError(er,describe);
    fprintf(PFH,"%s \n",describe);
  }
}
