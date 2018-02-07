#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_world_create(Fractal_Memory& mem,int level,vector <vector <int> >& SBoxes,
			  vector<vector<Point*>>& SPoints, bool buffer_groups)
  {
    static int _COUNTER=0;
    // int RANK=-1;
    // MPI_Comm_rank(Fractal_Memory::FRACTAL_UNIVERSE,&RANK);
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    mem.p_mess->IHranks.assign(FractalNodes,-1);
    clean_vector(mem.p_mess->Hranks);
    mem.p_mess->IAmAHypreNode=!SBoxes.empty();
    if(!buffer_groups)
      {
	if(mem.p_mess->IAmAHypreNode)
	  {
	    mem.p_mess->Hranks.push_back(FractalRank);
	    mem.p_mess->IHranks[FractalRank]=0;
	  }
	mem.p_mess->HypreNodes=mem.p_mess->Hranks.size();
	mem.p_mess->HypreGroupCreate(mem.p_mess->Hranks);
	return;
      }
    mem.Touchy.clear();
    int count=0;
    for(auto sp : SPoints)
      count+=sp.size();
    vector <int>counts(FractalNodes,0);
    mem.p_mess->How_Many_On_Nodes(count,counts);
    int countB=SBoxes.size();
    vector <int>countsB;
    if(mem.hypre_load_balance)
      {
	countsB.resize(FractalNodes,0);
	mem.p_mess->How_Many_On_Nodes(countB,countsB);
      }
    for(int FR : mem.TouchWhichBoxes)
      {
	if(counts[FR] > 0 && mem.p_mess->IAmAHypreNode)
	  {
	    vector <int>TBBox=mem.BBoxesLev[FR][level];
	    for(vector <int>& SB : SBoxes)
	      {
	    	if(!overlap_boxes(SB,TBBox))
	    	  continue;
	    	mem.Touchy.push_back(FR);
	    	break;
	      }
	  }
      }
    node_groups_struct(mem,counts);
    // mem.p_mess->freenodes.clear();
    if(mem.hypre_load_balance && !mem.p_mess->freenodes.empty())
      use_freenodes(mem,counts,countsB);
    mem.p_mess->HypreGroupCreate(mem.p_mess->Hranks);
    _COUNTER++;
  }
}
