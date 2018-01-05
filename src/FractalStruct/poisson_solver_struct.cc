#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver_struct(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    // time_t time_start=mem.p_mess->Clock();
    fractal.timing(-1,31);
    ofstream& FHT=mem.p_file->DUMPS;
    static int LEV=-1;
    static int _COUNTER=0;
    static int _COUNTERA=0;
    static vector<vector<int>> VOLA(2);
    static vector<vector<double>> FILLA(2);
    bool& IAmAHypreNode=mem.p_mess->IAmAHypreNode;
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    VOLA[0].resize(20);
    VOLA[1].resize(20);
    FILLA[0].resize(20);
    FILLA[1].resize(20);
    double time0,time1,time2,time3,time4,time5,time6,time7,time8;
    int spacing=Misc::pow(2,fractal.get_level_max()-level);
    for(int ni=0;ni<2;ni++)
      {
	mem.p_mess->IAmAHypreNode=mem.p_mess->count_on_node[2*level+ni];
	// FHT << " AA " << IAmAHypreNode << " " << level << " " << ni << "\n";
	// FHT.flush();
	bool doit=false;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	   // FHT << " COUNTS " << ni << " " << FR << " ";
	   // FHT << mem.p_mess->counts_on_nodes[ni+2*level][FR] << "\n";
	    if(!mem.p_mess->counts_on_nodes[ni+2*level][FR])
	      continue;
	    doit=true;
	    break;
	  }
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
	// FHT << " BB " << IAmAHypreNode << " " << level << " " << ni << "\n";
	// FHT.flush();
	time1=mem.p_mess->Clock();
	if(buffer)
	  hypre_points_clean(mem,level,hypre_points);
	// FHT << " CC " << IAmAHypreNode << " " << level << " " << ni << "\n";
	// FHT.flush();
	int VOLMIN=-100;
	double FILLFACTOR=100.7;
	int HPMIN=1000;
	if(!hypre_points.empty())
	  {
	    if(_COUNTERA % 10 == 0)
	      {
		hypre_best_boxes(mem,hypre_points,spacing,VOLMIN,FILLFACTOR);
		// FHT << " DD " << IAmAHypreNode << " " << level << " " << ni << "\n";
		// FHT.flush();
		VOLA[ni][level]=VOLMIN;
		FILLA[ni][level]=FILLFACTOR;
	      }
	    if(hypre_points.size() >= HPMIN)
	      hypre_points_boxes(mem,hypre_points,spacing,VOLA[ni][level],FILLA[ni][level],SBoxes,SPoints);
	    else
	      hypre_points_boxes(mem,hypre_points,spacing,VOLMIN,FILLFACTOR,SBoxes,SPoints);	      
	   // FHT << " EE " << IAmAHypreNode << " " << level << " " << ni << "\n";
	   // FHT.flush();
	    box_stats(mem,level,ni,SBoxes,SPoints);
	   // FHT << " FF " << IAmAHypreNode << " " << level << " " << ni << "\n";
	   // FHT.flush();
	  }
	time2=mem.p_mess->Clock();
	hypre_points.clear();
	hypre_world_create(mem,level,SBoxes,SPoints,buffer);
	// FHT << " GG " << IAmAHypreNode << " " << level << " " << ni << "\n";
	// FHT.flush();
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
	   // FHT << " HH " << IAmAHypreNode << " " << level << " " << ni << "\n";
	   // FHT.flush();
	    time5=mem.p_mess->Clock();
	    if(buffer)
	      add_buffer_values(mem,level,SBoxes,SPoints);
	   // FHT << " II " << IAmAHypreNode << " " << level << " " << ni << "\n";
	   // FHT.flush();
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
	// FHT << " JJ " << IAmAHypreNode << " " << level << " " << ni << "\n";
	// FHT.flush();
      }
    if(level < LEV)
      _COUNTERA++;
    _COUNTER++;
    LEV=level;
   // FHT << " Finished Poisson level " << level << "\n";
   // FHT.flush();
    fractal.timing(1,31);
  }
}
