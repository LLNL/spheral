#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_world_create(Fractal_Memory& mem,int level,vector <vector <int> >& SBoxes,
			  bool buffer_groups)
  {
//     int RANK=-1;
//     MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
//     bool Ranky = RANK == 21;
//     Ranky=true;
//     if(Ranky)
//       cerr << " Enter CREATE A " << RANK << " " << buffer_groups << endl;
    mem.p_mess->Full_Stop_Do_Not_Argue();
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    mem.p_mess->IHranks.assign(FractalNodes,-1);
    mem.p_mess->Hranks.clear();
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
    for(int FR : mem.TouchWhichBoxes)
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
    mem.p_mess->Full_Stop_Do_Not_Argue();
    node_groups_struct(mem);
    mem.p_mess->Full_Stop_Do_Not_Argue();
    mem.p_mess->HypreGroupCreate(mem.p_mess->Hranks);
    mem.p_mess->Full_Stop_Do_Not_Argue();
  }
}
