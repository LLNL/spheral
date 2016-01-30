#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_world_create(Fractal_Memory& mem,int level,vector <vector <int> >& SBoxes,
			 bool buffer_groups)
  {
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    mem.p_mess->IHranks.assign(FractalNodes,-1);
    mem.p_mess->Hranks.clear();
    mem.p_mess->IHranks.resize(FractalNodes);
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
    node_groups_struct(mem);
    mem.p_mess->HypreGroupCreate(mem.p_mess->Hranks);
  }
}
