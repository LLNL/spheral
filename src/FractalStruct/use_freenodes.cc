#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void use_freenodes(Fractal_Memory& mem,vector<int>& countsP,vector<int>& countsB)
  {
    if(mem.p_mess->freenodes.empty())
      return;
    ofstream& FHT=mem.p_file->DUMPS;
    const int atLEASTp=10000;
    const int maxLOADp=max(mem.hypre_max_node_load,atLEASTp);
    const int FractalRank=mem.p_mess->FractalRank;
    vector<int>sumPoints;
    vector<int>totalNodes;
    multimap<int,int>averagePoints;
    int number=0;
    mem.p_mess->mynumber=-1;
    for(auto nl : mem.p_mess->node_lists)
      {
	sumPoints.push_back(0);
	for(int FR : nl)
	  {
	    sumPoints.back()+=countsP[FR];
	    if(FR == FractalRank)
	      mem.p_mess->mynumber=number;
	  }
	totalNodes.push_back(nl.size());
	int aver=sumPoints.back()/totalNodes.back();
	if(aver > maxLOADp)
	  averagePoints.insert(make_pair(aver,number));
	FHT << " FREENODES A" << " " << mem.steps << " " << mem.level << " " << FractalRank << " " << number << " " << aver << " " << maxLOADp << "\n";
	number++;
      }
    FHT << " FREENODES B" << " " << mem.steps << " " << mem.level << " " << FractalRank << " " << averagePoints.size() << " " << mem.p_mess->freenodes.size() << "\n";
    // FHT.flush();
    if(averagePoints.empty())
      return;
    for(int FN : mem.p_mess->freenodes)
      {
	if(averagePoints.empty())
	  break;
	auto it=--averagePoints.end();
	int aver=it->first;
	if(aver <= maxLOADp)
	  break;
	int number=it->second;
	mem.p_mess->node_lists[number].push_back(FN);
	if(FN == FractalRank)
	  mem.p_mess->mynumber=number;
	averagePoints.erase(it);
	aver=(aver*(mem.p_mess->node_lists[number].size()-1))/mem.p_mess->node_lists[number].size();
	if(aver > maxLOADp)
	  averagePoints.insert(make_pair(aver,number));
	FHT << " FREENODES C" << mem.steps << " " << mem.level << " " << FractalRank << " " << FN << " " << number << " " << aver << " " << maxLOADp << "\n";
      }
    // FHT.flush();
    mem.p_mess->IHranks.clear();
    mem.p_mess->IHranks.resize(mem.p_mess->FractalNodes,-1);
    if(mem.p_mess->mynumber >= 0)
      {
	mem.p_mess->Hranks=mem.p_mess->node_lists[mem.p_mess->mynumber];
	int countit=0;
	for(auto HR : mem.p_mess->Hranks)
	  mem.p_mess->IHranks[HR]=countit++;
	mem.p_mess->IAmAHypreNode=true;
	mem.p_mess->HypreNodes=mem.p_mess->Hranks.size();
      }
    else
      {
	mem.p_mess->Hranks.clear();
	mem.p_mess->IAmAHypreNode=false;
	mem.p_mess->HypreNodes=0;
      }
  }
}
