#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_world_create(Fractal_Memory& mem,int level,vector <vector <int> >& SBoxes,
			  bool buffer_groups)
  {
    int RANK=-1;
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
    bool Ranky = RANK == 21;
    Ranky=true;
    if(Ranky)
      cerr << " CREATE A " << RANK << endl;
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
    int TWB=mem.TouchWhichBoxes.size();
    mem.Touchy.clear();
    for(int TW=0;TW<TWB;TW++)
      {
	int FR=mem.TouchWhichBoxes[TW];
	vector <int>TBBox=mem.BBoxesLev[FR][level];
	for(int B=0;B<SBoxes.size();B++)
	  {
	    if(!overlap_boxes(SBoxes[B],TBBox))
	      continue;
	    mem.Touchy.push_back(FR);
	    break;
	  }
      }
    if(Ranky)
      cerr << " Enter Node Groups 0 " << RANK << endl;
    mem.p_mess->Full_Stop_Do_Not_Argue();
    node_groups_struct(mem);
    if(Ranky)
      cerr << " Groups Create A " << RANK << endl;
    mem.p_mess->Full_Stop_Do_Not_Argue();
    mem.p_mess->HypreGroupCreate(mem.p_mess->Hranks);
    if(Ranky)
      cerr << " Groups Create B " << RANK << endl;
    mem.p_mess->Full_Stop_Do_Not_Argue();
  }
}
