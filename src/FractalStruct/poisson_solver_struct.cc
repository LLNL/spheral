#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver_struct(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    fractal.timing(-1,31);
    ofstream& FHT=mem.p_file->DUMPS;
    static int LEV=-1;
    static int _COUNTER=0;
    static int _COUNTERA=0;
    static vector<vector<int>> VOLA(2);
    static vector<vector<double>> FILLA(2);
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    VOLA[0].resize(20);
    VOLA[1].resize(20);
    FILLA[0].resize(20);
    FILLA[1].resize(20);
    //
    // vector<vector<int>> WorldRanks;
    // vector<vector<int>> LocalGroups;
    // vector<vector<int>> FreeNodes;
    // vector<bool> IAmIn;
    // double tsa=-mem.p_mess->Clock();
    // super_groups(mem,mem.all_groups[level],level,WorldRanks,LocalGroups,
    // 		 FreeNodes,IAmIn);
    // tsa+=mem.p_mess->Clock();
    // FHT << " SUPER timing " << tsa << " " << level << "\n";
    // WorldRanks.clear();
    // LocalGroups.clear();
    // FreeNodes.clear();
    // clean_vector(IAmIn);
    //
    double time0,time1,time2,time3,time4,time5,time6,time7,time8;
    int spacing=Misc::pow(2,fractal.get_level_max()-level);
    for(int ni=0;ni<2;ni++)
      {
	cerr.flush();
	mem.p_mess->IAmAHypreNode=mem.p_mess->count_on_node[2*level+ni];

	bool doit=false;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    if(!mem.p_mess->counts_on_nodes[ni+2*level][FR])
	      continue;
	    doit=true;
	    break;
	  }
	// cerr << " POISSON SOLVER ENTER " << mem.steps << " " << FractalRank << " " << level << " " << ni << " " << mem.p_mess->IAmAHypreNode << " " << doit << "\n";
	if(!doit)
	  continue;
	bool buffer=ni > 0;
	vector <vector <Point*> >hypre_points;
	vector <vector <Point*> >SPoints;
	vector <vector <int> >SBoxes;
	time0=mem.p_mess->Clock();
	bool single=buffer;
	single=false;
	hypre_points_struct(single,mem,mem.all_groups[level],hypre_points,buffer,level);
	time1=mem.p_mess->Clock();
	if(buffer)
	  hypre_points_clean(mem,level,hypre_points);
	int VOLMIN=-100;
	double FILLFACTOR=100.7;
	if(!hypre_points.empty())
	  {
	    if(_COUNTERA % 10 == 0)
	      {
		hypre_best_boxes(mem,hypre_points,spacing,VOLMIN,FILLFACTOR);
		VOLA[ni][level]=VOLMIN;
		FILLA[ni][level]=FILLFACTOR;
	      }
	    // FHT << " PUZZLED " << _COUNTERA << " " << ni << " " << level << " " << VOLA[ni][level] ;
	    // FHT << " " << FILLA[ni][level] << "\n";
	    hypre_points_boxes(mem,hypre_points,spacing,VOLA[ni][level],FILLA[ni][level],SBoxes,SPoints);
	    box_stats(mem,level,ni,SBoxes,SPoints);
	  }
	time2=mem.p_mess->Clock();
	hypre_points.clear();
	hypre_world_create(mem,level,SBoxes,SPoints,buffer);
	time3=mem.p_mess->Clock();
	double tt=-mem.p_mess->Clock();
	tt+=mem.p_mess->Clock();
	time4=time3;
	time5=time3;
	time6=time3;
	if(mem.p_mess->IAmAHypreNode)
	  {
	    time4=mem.p_mess->Clock();
	    hypre_solve_struct(ni > 0,mem,level,SBoxes,SPoints);
	    time5=mem.p_mess->Clock();
	    if(buffer)
	      add_buffer_values(mem,level,SBoxes,SPoints);
	    time6=mem.p_mess->Clock();
	  }
	time7=mem.p_mess->Clock();
	mem.p_mess->HypreGroupFree();
	time8=mem.p_mess->Clock();
	if(mem.p_mess->IAmAHypreNode)
	  {
	    FHT << " HYPRE RES B " <<  FractalRank << " " << ni << " " << level << " " << _COUNTERA << " ";
	    FHT << " " << time1-time0 << " " << time2-time1 << " " << time3-time2 << " " << time5-time4 << " " << time6-time5 <<  " " << time8-time7 << " " << tt << "\n";
	  }
	_COUNTER++;
	// cerr << " POISSON SOLVER EXIT " << mem.steps << " " << FractalRank << " " << level << " " << ni << "\n";
	// cerr.flush();
      }
    if(level < LEV)
      _COUNTERA++;
    _COUNTER++;
    LEV=level;
    fractal.timing(1,31);
  }
}
