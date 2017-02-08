#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool hypre_struct_load_balance(Fractal_Memory& mem,
				 vector<vector<int>>& SBoxes,
				 vector<vector<Point*>>& SPoints,
				 vector<int>& HRout)
  {
    static const int maxtries=4;
    const int maxLOAD=mem.hypre_max_node_load;
    const int extra=500;
    const double spread=0.1;
    ofstream& FHT=mem.p_file->DUMPS;
    int FractalRank=mem.p_mess->FractalRank;
    int HypreRank=mem.p_mess->HypreRank;
    int HypreNodes=mem.p_mess->HypreNodes;
    HRout.assign(SBoxes.size(),HypreRank);
    if(!mem.hypre_load_balance)
      return false;
    int count_points=0;
    for(auto SP : SPoints)
      count_points+=SP.size();
    vector<int>how_many;
    how_many.push_back(SPoints.size());
    how_many.push_back(count_points);
    vector<int>points_on_nodes(HypreNodes*2);
    mem.p_mess->my_AllgatherI(mem.p_mess->HypreWorld,how_many,points_on_nodes,2);
    vector<int>Boxes(HypreNodes);
    vector<int>Points(HypreNodes);
    int ni=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	Boxes[HR]=points_on_nodes[ni++];
	Points[HR]=points_on_nodes[ni++];
      }
    int total_boxes=accumulate(Boxes.begin(),Boxes.end(),0);
    int total_points=accumulate(Points.begin(),Points.end(),0);
    int average_boxes=total_boxes/HypreNodes;
    int average_points=total_points/HypreNodes;
    FHT << " LOADY " << total_boxes << " " << total_points << " " << average_boxes << " " << average_points << endl;
    multimap<int,deque<int>>NodesA;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	deque<int>VHR{{HR}};
	NodesA.insert(pair<int,deque<int>>(Points[HR],VHR));
      }
    for(auto M : NodesA)
      FHT << " CEND " << FractalRank << " " << HypreRank << " " << mem.steps << " " << mem.level << " " << M.first << " " << M.second.front() << " " << maxLOAD << endl;
    if((--NodesA.end())->first <= maxLOAD)
      return false;
    int enough_spam=false;
    int tries=0;
    while((--NodesA.end())->first > maxLOAD && tries < maxtries && !enough_spam)
      {
	multimap<int,deque<int>>NodesB;
	auto pstart=NodesA.begin();
	auto pend=--NodesA.end();
	int dd=distance(pstart,pend);
	while(dd > 0)
	  {
	    FHT << " HERE 0 " << distance(pstart,pend) << " " << dd << " " << FractalRank << endl;
	    auto Va=pstart->second;
	    auto Vb=pend->second;
	    int aver=(Va.size()*pstart->first+Vb.size()*pend->first)/
	      (Va.size()+Vb.size());	  
	    FHT << " HERE A " << dd << " " << pstart->first << " " << pend->first << " " << aver << " " << " " << Va.size() << " " << Vb.size() << " " << FractalRank << endl;
	    for(auto V : Vb)
	      Va.push_back(V);
	    NodesB.insert(make_pair(aver,Va));
	    pstart++;
	    pend--;
	    dd-=2;
	  }
	if(dd == 0)
	  NodesB.insert(*pstart);
	NodesA=NodesB;
	FHT << " HERE C " << tries << " " << NodesA.size() << endl;
	tries++;
	enough_spam=(double)(NodesB.begin()->first-(--NodesB.end())->first)/
	  (double)(NodesB.begin()->first+(--NodesB.end())->first) < spread;
      }
    bool go_home=true;
    deque<int>mynodes;
    int aver=0;
    for(auto No : NodesA)
      {
	mynodes=No.second;
	if(find(mynodes.begin(),mynodes.end(),HypreRank) != mynodes.end())
	  {
	    for(auto HR : mynodes)
	      {
		FHT << " Here DD " << HR << " " << Points[HR] << endl;
		aver+=Points[HR];
		go_home=go_home && Points[HR] <= maxLOAD;
	      }
	    break;
	  }
      }
    aver/=mynodes.size();
    if(go_home)
      return true;
    NodesA.clear();
    int averup=(aver*102)/100;
    int averdown=(aver*98)/100;
    vector<int>can(HypreNodes,0);
    vector<int>need(HypreNodes,0);
    sort(mynodes.begin(),mynodes.end());
    for(auto HR : mynodes)
      {
	can[HR]=max(Points[HR]-averup,0);
	need[HR]=max(averdown-Points[HR],0);
	FHT << " CAN NEED A " << HR << " " << can[HR] << " " << need[HR] << endl;
      }
    vector<int>sendto(HypreNodes,0);
    for(int Hsend : mynodes)
      {
	if(Hsend > HypreRank)
	  break;
	if(can[Hsend] <= 0)
	  continue;
	for (int Hrec : mynodes)
	  {
	    if(need[Hrec] <= 0 || Hrec == Hsend)
	      continue;
	    int send=min(can[Hsend],need[Hrec]);
	    FHT << " CAN NEED B " << Hsend << " " << Hrec << " " << can[Hsend] << " " << need[Hrec] << " " << send << endl;
	    can[Hsend]-=send;
	    need[Hrec]-=send;
	    if(Hsend == HypreRank)
	      sendto[Hrec]=send;
	  }
      }
    FHT << " MYNODES A " << mynodes.size();
    for(auto pN=mynodes.begin();pN!=mynodes.end();pN++)
      if(*pN == HypreRank)
	{
	  mynodes.erase(pN);
	  break;
	}
    FHT << " " << mynodes.size() << endl;
    sendto[HypreRank]=-1;
    int number=0;
    multimap <vector<Point*>,int,vector_comp_down>MSP;
    for(auto SP : SPoints)
      MSP.insert(make_pair(SP,number++));
    bool fewer=false;
    for(auto M : MSP)
      {
	int has=M.first.size();
	FHT << " MYNODES B " << MSP.size() << " " << has;
	for(auto Hrec : mynodes)
	  {
	    FHT << " " << Hrec << " " << sendto[Hrec];
	    if(sendto[Hrec] > 0)
	      {
		if(has <= sendto[Hrec]+extra)
		  {
		    HRout[M.second]=Hrec;
		    FHT << " " << M.second << " " << Hrec;;
		    sendto[Hrec]-=has;
		    break;
		  }
	      }
	    else
	      fewer=true;
	  }
	FHT << endl;
	if(!fewer)
	  continue;
	auto pN=mynodes.begin();
	while(pN != mynodes.end())
	  {
	    if(sendto[*pN] <= 0)
	      pN=mynodes.erase(pN);
	    else
	      pN++;
	  }
      }
    int spam=0;
    for(auto SP : SPoints)
      {
	if(HRout[spam] != HypreRank)
	  FHT << " SENDIT A " << HypreRank << " " << HRout[spam] << " " << SP.size() << endl;
	spam++;
      }
    return true;
  }
}
