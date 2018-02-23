#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_solve_struct(bool buffer,Fractal_Memory& mem,int level,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    static int _COUNTER=0;
    static vector <double> Hypre_sum_time(mem.p_fractal->get_level_max()+1,0.0);
    ofstream& FHT=mem.p_file->DUMPS;
    int spacing=Misc::pow(2,mem.p_fractal->get_level_max()-level);
    int FractalRank=mem.p_mess->FractalRank;
    int HypreRank=mem.p_mess->HypreRank;
    int HypreNodes=mem.p_mess->HypreNodes;
    vector <int> counts_in(HypreNodes,0);
    vector <int> counts_out(HypreNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(HypreNodes);
    vector < vector <double> > dataR_out(HypreNodes);

    vector<int>HRout(SBoxes.size(),HypreRank);
    int balance=false;
    if(buffer && mem.hypre_load_balance)
      balance=hypre_struct_load_balance(mem,SBoxes,SPoints,HRout);
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   Amatrix;
    HYPRE_StructVector   rho;
    HYPRE_StructVector   pot;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;
    vector<int>pos(3);
    vector <int>Box=mem.FRBoxesLev[level];
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
    int Bstay=0;
    int Bout=0;
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
	if(HRout[B]==HypreRank)
	  {
	    HYPRE_StructGridSetExtents(grid,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()));
	    sumVOL+=VOL[B];
	    Bstay++;
	  }
	else
	  {
	    int HR=HRout[B];
	    for(int ni=0;ni<3;ni++)
	      dataI_out[HR].push_back(lowerBOX[B][ni]);
	    for(int ni=0;ni<3;ni++)
	      dataI_out[HR].push_back(upperBOX[B][ni]);
	    counts_out[HR]++;
	    Bout++;
	  }
	B++;
      }
    assert(Bstay+Bout==B);
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=6;
    int doubles=0;
    double timeM0=-mem.p_mess->Clock();
    if(balance)
      mem.p_mess->Send_Data_Some_How(38,mem.p_mess->HypreWorld,
				     counts_out,counts_in,integers,doubles,
				     dataI_out,dataI_in,how_manyI,
				     dataR_out,dataR_in,how_manyR);
    clean_vector(counts_out);
    dataI_out.clear();
    dataR_out.clear();
    int c6=0;
    int Bget=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	for(int c=0;c<counts_in[HR];c++)
	  {
	    HRout.push_back(-HR-1);
	    lowerBOX.resize(lowerBOX.size()+1);
	    upperBOX.resize(upperBOX.size()+1);
	    VOL.push_back(1);
	    for(int ni : {0,1,2})
	      lowerBOX.back().push_back(dataI_in[c6++]);
	    for(int ni : {0,1,2})
	      {
		upperBOX.back().push_back(dataI_in[c6++]);
		VOL.back()*=upperBOX.back()[ni]-lowerBOX.back()[ni]+1;
		assert(VOL.back() > 0);
	      }
	    HYPRE_StructGridSetExtents(grid,&(*lowerBOX.back().begin()),&(*upperBOX.back().begin()));
	    sumVOL+=VOL.back();
	    Bget++;
	  }
      }
    clean_vector(dataI_in);
    clean_vector(dataR_in);
    clean_vector(counts_in);
    timeM0+=mem.p_mess->Clock();
    HYPRE_StructGridAssemble(grid);
    long int SVT=mem.p_mess->How_Many_In_Solver(sumVOL);
    long int SBT=mem.p_mess->How_Many_In_Solver(Bstay+Bget);
    double time1=mem.p_mess->Clock();
    
    HYPRE_StructStencilCreate(3,7,&stencil);
    int offsets[7][3] = {{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    for (int entry = 0; entry < 7; entry++)
      HYPRE_StructStencilSetElement(stencil,entry,offsets[entry]);

    counts_in.resize(HypreNodes,0);
    counts_out.resize(HypreNodes,0);
    dataI_out.resize(HypreNodes);
    dataR_out.resize(HypreNodes);
    
    HYPRE_StructMatrixCreate(mem.p_mess->HypreWorld,grid,stencil, &Amatrix);
    HYPRE_StructMatrixInitialize(Amatrix);
    HYPRE_StructVectorCreate(mem.p_mess->HypreWorld, grid, &rho);
    HYPRE_StructVectorCreate(mem.p_mess->HypreWorld, grid, &pot);
    HYPRE_StructVectorInitialize(rho);
    HYPRE_StructVectorInitialize(pot);
    int stencil_indices[7] = {0,1,2,3,4,5,6};
    B=0;
    for(auto& SP : SPoints)
      {
	vector <double>values;
	vector <double>dens_values;
	vector <double>pot_values;
	for(auto p : SP)
	  {
	    if(p == 0)
	      {
		values.push_back(1.0);
		for(int ni=0;ni<6;ni++)
		  values.push_back(0.0);
		dens_values.push_back(1.0);
		pot_values.push_back(1.0);
	      }
	    else
	      {
		values.push_back(-6.0);
		for (int ni = 0; ni < 6; ni++)
		  {
		    Point* p1=p->get_point_ud_0(ni);
		    if(p1->get_inside() && !p1->get_trouble())
		      values.push_back(1.0);
		    else
		      values.push_back(0.0);
		  }
		double density=p->get_density_point()*g_c;
		for(int ni=0;ni<6;ni++)
		  {
		    Point* p1=p->get_point_ud_0(ni);
		    if(!p1->get_inside() || p1->get_trouble())
		      density-=p1->get_potential_point();
		  }
		dens_values.push_back(density);
		pot_values.push_back(p->get_potential_point());
	      }
	  }
	if(HRout[B] == HypreRank)
	  {
	    HYPRE_StructMatrixAddToBoxValues(Amatrix,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),7,
					     stencil_indices,&(*values.begin()));
	    HYPRE_StructVectorAddToBoxValues(rho,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),&(*dens_values.begin()));
	    HYPRE_StructVectorAddToBoxValues(pot,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),&(*pot_values.begin()));
	  }
	else
	  {
	    int HR=HRout[B];
	    int knights=dens_values.size();
	    int ni7=0;
	    for(int ni=0;ni<knights;ni++)
	      {
		dataR_out[HR].push_back(dens_values[ni]);
		dataR_out[HR].push_back(pot_values[ni]);
		int NV=0;
		for(int ni=6;ni>0;ni--)
		  {
		    int nn=(int)(values[ni7+ni]+0.01);
		    NV=2*NV+nn;
		  }
		NV*=2;
		int n0=(int)(values[ni7]+0.01);
		if(n0 >= 0)
		  NV++;
		dataI_out[HR].push_back(NV);
		ni7+=7;
		counts_out[HR]++;
	      }
	  }
	B++;
      }
    how_manyI=-1;
    how_manyR=-1;
    integers=1;
    doubles=2;
    double timeM1=-mem.p_mess->Clock();
    if(balance)
      mem.p_mess->Send_Data_Some_How(48,mem.p_mess->HypreWorld,
				     counts_out,counts_in,integers,doubles,
				     dataI_out,dataI_in,how_manyI,
				     dataR_out,dataR_in,how_manyR);
    clean_vector(counts_out);
    dataI_out.clear();
    dataR_out.clear();
    int ni=0;
    int nd=0;
    B=Bstay+Bout;
    for(int Bg=0;Bg<Bget;Bg++)
      {
	int vol=1;
	for(int knights : {0,1,2})
	  vol*=upperBOX[B][knights]-lowerBOX[B][knights]+1;
	assert(vol == VOL[B]);
	vector <double>values;
	vector <double>dens_values;
	vector <double>pot_values;
	for(int knights=0;knights<VOL[B];knights++)
	  {
	    dens_values.push_back(dataR_in[nd++]);
	    pot_values.push_back(dataR_in[nd++]);
	    int NV=dataI_in[ni++];
	    if(NV % 2 == 1)
	      values.push_back(1.0);
	    else
	      values.push_back(-6.0);
	    for(int V=0;V<6;V++)
	      {
		NV/=2;
		values.push_back((double)(NV % 2));
	      }
	  }
	HYPRE_StructMatrixAddToBoxValues(Amatrix,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),7,
					 stencil_indices,&(*values.begin()));
	HYPRE_StructVectorAddToBoxValues(rho,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),&(*dens_values.begin()));
	HYPRE_StructVectorAddToBoxValues(pot,&(*lowerBOX[B].begin()),&(*upperBOX[B].begin()),&(*pot_values.begin()));
	B++;
      }
    clean_vector(dataI_in);
    clean_vector(dataR_in);
    clean_vector(counts_in);
    timeM1+=mem.p_mess->Clock();
    HYPRE_StructMatrixAssemble(Amatrix);
    HYPRE_StructVectorAssemble(rho);
    HYPRE_StructVectorAssemble(pot);
    double time2=mem.p_mess->Clock();
    
    HYPRE_StructPCGCreate(mem.p_mess->HypreWorld, &solver);
    HYPRE_StructPCGSetMaxIter(solver, mem.maxits );
    HYPRE_StructPCGSetTol(solver, mem.HTOL );
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
    FHT << " SOLVED A " << level << " " << _COUNTER << " " << FractalRank << " " << HypreRank << " " << num_iterations << " " << final_res_norm << "\n";
    HYPRE_StructPCGDestroy(solver);
    HYPRE_StructPFMGDestroy(precond);
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructMatrixDestroy(Amatrix);
    HYPRE_StructVectorDestroy(rho);
    counts_in.resize(HypreNodes,0);
    counts_out.resize(HypreNodes,0);
    dataI_out.resize(HypreNodes);
    dataR_out.resize(HypreNodes);
    B=0;
    for(int BOX=0;BOX<lowerBOX.size();BOX++)
      {
	vector <double>pot_values;
	pot_values.resize(VOL[BOX]);
	if(HRout[BOX] == HypreRank)
	  {
	    HYPRE_StructVectorGetBoxValues(pot,&(*lowerBOX[BOX].begin()),&(*upperBOX[BOX].begin()),&(*pot_values.begin()));
	    int ni=0;
	    for(Point* &p : SPoints[BOX])
	      {
		if(p != 0 && !p->get_trouble())
		  p->set_potential_point(pot_values[ni]);
		ni++;
	      }
	  }
	else if(HRout[BOX] < 0)
	  {
	    HYPRE_StructVectorGetBoxValues(pot,&(*lowerBOX[BOX].begin()),&(*upperBOX[BOX].begin()),&(*pot_values.begin()));
	    int HR=-HRout[BOX]-1;
	    for(int ni=0;ni<VOL[BOX];ni++)
	      {
		dataR_out[HR].push_back(pot_values[ni]);
		counts_out[HR]++;
	      }
	  }
	B++;
      }
    HYPRE_StructVectorDestroy(pot);
    how_manyI=-1;
    how_manyR=-1;
    integers=0;
    doubles=1;
    double timeM2=-mem.p_mess->Clock();
    if(balance)
      mem.p_mess->Send_Data_Some_How(58,mem.p_mess->HypreWorld,
				     counts_out,counts_in,integers,doubles,
				     dataI_out,dataI_in,how_manyI,
				     dataR_out,dataR_in,how_manyR);
    clean_vector(counts_out);
    dataI_out.clear();
    dataR_out.clear();
    int HR=-1;
    ni=0;
    for(auto &SP : SPoints)
      {
	HR++;
	if(HRout[HR] < 0 || HRout[HR] == HypreRank)
	  continue;
	for(auto &p : SP)
	  {
	    double pot=dataR_in[ni++];
	    if(p != 0 && !p->get_trouble())
	      p->set_potential_point(pot);
	  }
      }
    timeM2+=mem.p_mess->Clock();
    double time6=mem.p_mess->Clock();
    _COUNTER++;
    Hypre_sum_time[level]+=time6-time0;
    double timeM=timeM0+timeM1+timeM2;
    FHT << " Hypre Total " << FractalRank << " " << time6-time0 << " " << Hypre_sum_time[level] << " L " << level << " B " << buffer << " " << sumVOL << " " << SVT << " " << SBoxes.size() << " " << SBT << " S " << mem.steps << "\n";
    FHT << " Hypre Grid Assemble " << "\t" << time1-time0 << "\n";
    FHT << " Hypre Data Assemble " << "\t" << time2-time1 << "\n";
    FHT << " Hypre Solve Assemble " << "\t" << time3-time2 << "\n";
    FHT << " Hypre Solve Setup " << "\t" << time4-time3 << "\n";
    FHT << " Hypre Solver Solve " << "\t" << time5-time4 << "\n";
    FHT << " Hypre Data Dump " << "\t" << time6-time5 << "\n";
    FHT << " Hypre MPI Stuff " << "\t" << timeM  << "\t" << timeM0 << "\t" << timeM1 << "\t" << timeM2 << "\n";
    if(buffer)
      {
	int Pstay=0;
	int Pout=0;
	int Pget=0;
	for(int ni=0;ni<HRout.size();ni++)
	  {
	    if(HRout[ni] == HypreRank)
	      Pstay+=VOL[ni];
	    else if(HRout[ni] < 0)
	      Pget+=VOL[ni];
	    else
	      Pout+=VOL[ni];
	  }
	FHT << " NODE COUNTS " << mem.level << " " << mem.steps << " " << mem.p_mess->mynumber << " " << FractalRank << " " << HypreRank << " " << Bstay+Bout << " " << Bstay+Bget << " " << Bstay << " " << Bout << " " << Bget <<"\n";
	FHT << " POINT COUNTS " << mem.level << " " << mem.steps << " " << mem.p_mess->mynumber << " " << FractalRank << " " << HypreRank << " " << Pstay+Pout << " " << Pstay+Pget << " " << Pstay << " " << Pout << " " << Pget <<"\n";
      }
  }
}
