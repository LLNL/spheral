#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_solve_struct(Fractal_Memory& mem,int level,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    static int _COUNTER=0;
    static vector <double> Hypre_sum_time(mem.p_fractal->get_level_max()+1,0.0);
    ofstream& FHT=mem.p_file->DUMPS;
    int spacing=Misc::pow(2,mem.p_fractal->get_level_max()-level);
    int FractalRank=mem.p_mess->FractalRank;
    int HypreRank=mem.p_mess->HypreRank;
    // int HypreNodes=mem.p_mess->HypreNodes;
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   Amatrix;
    HYPRE_StructVector   rho;
    HYPRE_StructVector   pot;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;
    vector<int>pos(3);
    vector <int>Box=mem.FRBoxesLev[level];
    FHT << " HYP BOX " << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3] << " " << Box[4] << " " << Box[5] << "\n";
    double pi = 4.0*atan(1.0);
    int length=mem.p_fractal->get_grid_length();
    double g_c=4.0*pi/static_cast<double>(length*length)*pow(4.0,-level);
    vector < vector <int> > lowerBOX(SBoxes.size());
    vector < vector <int> > upperBOX(SBoxes.size());
    vector <int> VOL(SBoxes.size());
    double time0=mem.p_mess->Clock();
    HYPRE_StructGridCreate(mem.p_mess->HypreWorld,3,&grid);
    vector <int>pers(3,0);
    HYPRE_StructGridSetPeriodic(grid,&(*pers.begin()));
    int sumVOL=0;
    int B=0;
    for(vector <int>& SB : SBoxes)
      {
	int ni=0;
	VOL[B]=1;
	for(int ni2=0;ni2<6;ni2+=2)
	  {
	    lowerBOX[B].push_back(SB[ni2]/spacing);
	    upperBOX[B].push_back(SB[ni2+1]/spacing);
	    VOL[B]*=upperBOX[B][ni]-lowerBOX[B][ni]+1;
	    assert(VOL[B] > 0);
	    ni++;
	  }
	HYPRE_StructGridSetExtents(grid,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()));
	sumVOL+=VOL[B];
	B++;
      }
    HYPRE_StructGridAssemble(grid);
    long int SVT=mem.p_mess->How_Many_In_Solver(sumVOL);
    long int SBT=mem.p_mess->How_Many_In_Solver(SBoxes.size());
    double time1=mem.p_mess->Clock();
    HYPRE_StructStencilCreate(3,7,&stencil);
    int offsets[7][3] = {{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    for (int entry = 0; entry < 7; entry++)
      HYPRE_StructStencilSetElement(stencil,entry,offsets[entry]);
    HYPRE_StructMatrixCreate(mem.p_mess->HypreWorld,grid,stencil, &Amatrix);
    HYPRE_StructMatrixInitialize(Amatrix);
    int stencil_indices[7] = {0,1,2,3,4,5,6};
    B=0;
    for(vector <Point*>& SP : SPoints)
      {
	vector <double>values;
	auto p_itr=SP.begin();
	while(p_itr != SP.end())
	  {
	    Point* p=*p_itr;
	    if(p != 0)
	      {
		p->get_pos_point(pos);
		values.push_back(-6.0);
		for (int j = 0; j < 6; j++)
		  {
		    Point* p1=p->get_point_ud_0(j);
		    p1->get_pos_point(pos);
		    bool good=p1->get_inside();
		    if(good)
		      {
			good=vector_in_box(pos,Box);
			if(!good)
			  {
			    for(auto FR : mem.Touchy)
			      {
				good=vector_in_box(pos,mem.BoxesLev[FR][level]);
				if(good)
				  break;
			      }
			  }
		      }
		    if(good)
		      values.push_back(1.0);
		    else
		      values.push_back(0.0);
		  }
	      }
	    else
	      {
		values.push_back(1.0);
		for(int ni=0;ni<6;ni++)
		  values.push_back(0.0);
	      }
	    p_itr++;
	  }
	HYPRE_StructMatrixAddToBoxValues(Amatrix,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),7,
					 stencil_indices,&(*values.begin()));
	values.clear();
	B++;
      }
    HYPRE_StructMatrixAssemble(Amatrix);
    HYPRE_StructVectorCreate(mem.p_mess->HypreWorld, grid, &rho);
    HYPRE_StructVectorCreate(mem.p_mess->HypreWorld, grid, &pot);
    
    HYPRE_StructVectorInitialize(rho);
    HYPRE_StructVectorInitialize(pot);
    B=0;
    for(vector <Point*>& SP : SPoints)
      {
	vector <double>dens_values;
	vector <double>pot_values;
	for(auto &p : SP)
	  {
	    if(p == 0)
	      {
		dens_values.push_back(1.0);
		pot_values.push_back(1.0);
	      }
	    else
	      {
		double density=p->get_density_point()*g_c;
		for(int ni=0;ni<6;ni++)
		  {
		    Point* p1=p->get_point_ud_0(ni);
		    bool good=p1->get_inside();
		    if(good)
		      {
			p1->get_pos_point(pos);
			good=vector_in_box(pos,Box);
			if(!good)
			  {
			    for(auto FR : mem.Touchy)
			      {
				good=vector_in_box(pos,mem.BoxesLev[FR][level]);
				if(good)
				  break;
			      }
			  }
		      }
		    if(!good)
		      density-=p1->get_potential_point();
		  }
		dens_values.push_back(density);
		pot_values.push_back(p->get_potential_point());
	      }
	  }
	HYPRE_StructVectorAddToBoxValues(rho,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),&(*dens_values.begin()));
	HYPRE_StructVectorAddToBoxValues(pot,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),&(*pot_values.begin()));
	dens_values.clear();
	pot_values.clear();
	B++;
      }
    HYPRE_StructVectorAssemble(rho);
    HYPRE_StructVectorAssemble(pot);
    double time2=mem.p_mess->Clock();
    
    HYPRE_StructPCGCreate(mem.p_mess->HypreWorld, &solver);
    HYPRE_StructPCGSetMaxIter(solver, 200 );
    HYPRE_StructPCGSetTol(solver, 1.0e-06 );
    HYPRE_StructPCGSetTwoNorm(solver, 1 );
    HYPRE_StructPCGSetRelChange(solver, 0 );
    HYPRE_StructPCGSetPrintLevel(solver, 1 );
  
    HYPRE_StructPFMGCreate(mem.p_mess->HypreWorld, &precond);
    HYPRE_StructPFMGSetMaxIter(precond, 1);
    HYPRE_StructPFMGSetTol(precond, 0.0);
    HYPRE_StructPFMGSetRAPType(precond, 1);
    HYPRE_StructPFMGSetSkipRelax(precond, 1);
    HYPRE_StructPFMGSetPrintLevel(precond, 0);
    HYPRE_StructPFMGSetLogging(precond, 0);
    HYPRE_StructPCGSetPrecond(solver,
			      HYPRE_StructPFMGSolve,
			      HYPRE_StructPFMGSetup,
			      precond);
    double time3=mem.p_mess->Clock();
    HYPRE_StructPCGSetup(solver,Amatrix,rho,pot);
    double time4=mem.p_mess->Clock();
    HYPRE_StructPCGSolve(solver,Amatrix,rho,pot);
    double time5=mem.p_mess->Clock();
    int num_iterations=-1;
    double final_res_norm=-1.0;
    HYPRE_StructPCGGetNumIterations(solver,&num_iterations );
    HYPRE_StructPCGGetFinalRelativeResidualNorm( solver, &final_res_norm );
    // if(mem.p_mess->IAmAHypreNode && HypreRank == 0)
    FHT << " SOLVED A " << level << " " << _COUNTER << " " << FractalRank << " " << HypreRank << " " << num_iterations << " " << final_res_norm << "\n";
    HYPRE_StructPCGDestroy(solver);
    HYPRE_StructPFMGDestroy(precond);
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructMatrixDestroy(Amatrix);
    HYPRE_StructVectorDestroy(rho);
    // cerr << " SOLVED B " << _COUNTER << " " << FractalRank << " " << HypreRank << "\n";
    B=0;
    for(vector <Point*>& SP : SPoints)
      {
	vector <double>pot_values;
	pot_values.resize(VOL[B]);
	HYPRE_StructVectorGetBoxValues(pot,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),&(*pot_values.begin()));
	int i=0;
	for(Point* &p : SP)
	  {
	    if(p != 0)
	      p->set_potential_point(pot_values[i]);
	    i++;
	  }
	B++;
      }
    HYPRE_StructVectorDestroy(pot);
    double time6=mem.p_mess->Clock();
    _COUNTER++;
    Hypre_sum_time[level]+=time6-time0;
    FHT << " Hypre Total " << FractalRank << " " << time6-time0 << " " << Hypre_sum_time[level] << " L" << level << " " << sumVOL << " " << SVT << " " << SBoxes.size() << " " << SBT << " S" << mem.steps << "\n";
    FHT << " Hypre Grid Assemble " << "\t" << time1-time0 << "\n";
    FHT << " Hypre Data Assemble " << "\t" << time2-time1 << "\n";
    FHT << " Hypre Solve Assemble " << "\t" << time3-time2 << "\n";
    FHT << " Hypre Solve Setup " << "\t" << time4-time3 << "\n";
    FHT << " Hypre Solver Solve " << "\t" << time5-time4 << "\n";
    FHT << " Hypre Data Dump " << "\t" << time6-time5 << "\n";
    // cerr << " SOLVED C " << _COUNTER << " " << FractalRank << " " << HypreRank << "\n";
  }
}
