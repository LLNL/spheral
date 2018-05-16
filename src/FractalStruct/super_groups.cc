#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void super_groups(Fractal_Memory& mem,vector <Group*>& groups,const int level,
		    vector<vector<int>>& WorldRanks,
		    vector<vector<int>>& LocalGroups,
		    vector<vector<int>>& FreeNodes,
		    vector<bool>& IAmIn)
  {
    ofstream& FHT=mem.p_file->DUMPS;
    clean_vector(mem.Touchy);
    WorldRanks.clear();
    LocalGroups.clear();
    FreeNodes.clear();
    clean_vector(IAmIn);
    clean_vector(mem.Touchy);
    WorldRanks.push_back(vector<int>());
    LocalGroups.push_back(vector<int>());
    FreeNodes.push_back(vector<int>());
    IAmIn.push_back(false);
    const int FractalNodes=mem.p_mess->FractalNodes;
    const int FractalRank=mem.p_mess->FractalRank;
    vector<int>LowTouch;
    vector<vector<int>>LowBox;
    vector <int>MyBox(mem.BoxesLev[FractalRank][level]);
    int counta=0;
    int countb=0;
    for(auto pg : groups)
      if(pg->get_buffer_group())
	countb=2;
      else
	counta=1;
    vector<int>counts(FractalNodes);
    mem.p_mess->How_Many_On_Nodes(counta+countb,counts);
    vector<bool>occupieda;
    vector<bool>occupiedb;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	occupieda.push_back(counts[FR] % 2 == 1);
	occupiedb.push_back((counts[FR]/2) % 2 == 1);
      }
    for(auto FR : mem.TouchWhichBoxes)
      {
	if(!occupiedb[FR] || !occupiedb[FractalRank])
	  continue;
	mem.Touchy.push_back(FR);
	if(FR >= FractalRank)
	  continue;
	LowTouch.push_back(FR);
	LowBox.push_back(mem.BoxesLev[FR][level]);
      }
    map<array<int,3>,int>edgies;
    vector <int> counts_in(FractalNodes);
    vector <int> dataI_in;
    vector <int> counts_out(FractalNodes,0);
    vector <double> dataR_in;
    vector <vector<int>> dataI_out(FractalNodes);
    vector <vector<double>> dataR_out(FractalNodes);
    int group_number=-1;
    vector<int>ggsize;
    for(auto pg : groups)
      {
	group_number++;
	if(!pg->get_buffer_group())
	  continue;
	int NP=0;
	for(auto p : pg->list_points)
	  {
	    if(!p->get_inside())
	      continue;
	    NP++;
	    if(!on_edge(p,MyBox))
	      continue;
	    edgies[p->get_pos_point_a()]=group_number;
	    for(int ni=0;ni<6;ni++)
	      {
		Point* pa=p->get_point_ud_0(ni,-27);
		if(vector_in_box(pa,MyBox))
		  continue;
		for(int FR : LowTouch)
		  {
		    if(!on_edge(pa,mem.BoxesLev[FR][level]))
		      continue;
		    array<int,3>pos(pa->get_pos_point_a());
		    dataI_out[FR].push_back(pos[0]);
		    dataI_out[FR].push_back(pos[1]);
		    dataI_out[FR].push_back(pos[2]);
		    dataI_out[FR].push_back(group_number);
		    counts_out[FR]++;
		    break;
		  }
	      }
	  }
	ggsize.push_back(group_number);
	ggsize.push_back(NP);
      }
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=4;
    int doubles=0;
    mem.p_mess->Send_Data_Some_How(12,mem.p_mess->FractalWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    clean_vector(counts_out);
    dataI_out.clear();
    dataR_out.clear();
    int c4=0;
    set<array<int,3>>GroupHug;
    for(int FR=0;FR<FractalNodes;FR++)
      for(int c=0;c<counts_in[FR];c++)
	{
	  auto it=edgies.find({{dataI_in[c4],dataI_in[c4+1],dataI_in[c4+2]}});
	  if(it != edgies.end())
	    GroupHug.insert({{it->second,FR,dataI_in[c4+3]}});
	  c4+=4;
	}
    clean_vector(dataI_in);
    clean_vector(dataR_in);
    clean_vector(counts_in);
    counts_out.resize(FractalNodes,0);
    dataI_out.resize(FractalNodes);
    int gg2=ggsize.size()/2;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	if(!occupiedb[FR])
	  continue;
	dataI_out[FR]=ggsize;
	counts_out[FR]=gg2;
      }
    how_manyI=-1;
    how_manyR=-1;
    integers=2;
    doubles=0;
    mem.p_mess->Send_Data_Some_How(22,mem.p_mess->FractalWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    clean_vector(counts_out);
    dataI_out.clear();
    dataR_out.clear();
    vector<vector<map<pair<int,int>,int>>>S2NGroups;
    vector<vector<set<int>>>S2Nodes;
    map<pair<int,int>,int>NodeGroups; // FR, Group#, Points in Group    
    int c2=0;
    for(int FR=0;FR<FractalNodes;FR++)
      for(int c=0;c<counts_in[FR];c++)
	{
	  NodeGroups[{FR,dataI_in[c2]}]=dataI_in[c2+1];
	  c2+=2;
	}
    clean_vector(dataI_in);
    clean_vector(dataR_in);
    clean_vector(counts_in);
    counts_out.resize(FractalNodes,0);
    dataI_out.resize(FractalNodes);
    dataR_out.resize(FractalNodes);
    how_manyI=-1;
    how_manyR=-1;
    integers=3;
    doubles=0;
    vector<int>DI;
    for(auto GH : GroupHug)
      {
	DI.push_back(GH[0]);
	DI.push_back(GH[1]);
	DI.push_back(GH[2]);
      }
    GroupHug.clear();
    const int co=DI.size()/3;
    for(int FR=0;FR<FractalNodes;FR++)
      if(occupiedb[FR])
	{
	  dataI_out[FR]=DI;
	  counts_out[FR]=co;
	}
    clean_vector(DI);
    mem.p_mess->Send_Data_Some_How(39,mem.p_mess->FractalWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    clean_vector(counts_out);
    dataI_out.clear();
    dataR_out.clear();
    clean_vector(dataR_in);
    vector<bool>NGfound(NodeGroups.size(),false);
    int c3=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    auto ita=NodeGroups.find({FR,dataI_in[c3]});
	    auto itb=NodeGroups.find({dataI_in[c3+1],dataI_in[c3+2]});
	    if(ita != NodeGroups.end() && itb != NodeGroups.end());
	    {
	      int nia=distance(NodeGroups.begin(),ita);
	      int nib=distance(NodeGroups.begin(),itb);
	      NGfound[nia]=true;
	      NGfound[nib]=true;
	    }
	    c3+=3;
	  }
      }
    auto pNGa=NodeGroups.begin();
    const auto pNGb=NodeGroups.end();
    int ni=0;
    while(pNGa != pNGb)
      {
	if(!NGfound[ni++])
	  pNGa=NodeGroups.erase(pNGa);
	else
	  pNGa++;
      }
    vector<int>list_pair_1;
    vector<int>list_pair_2;
    c3=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    auto it=NodeGroups.find({FR,dataI_in[c3]});
	    list_pair_1.push_back(distance(NodeGroups.begin(),it));
	    it=NodeGroups.find({dataI_in[c3+1],dataI_in[c3+2]});
	    list_pair_2.push_back(distance(NodeGroups.begin(),it));
	    c3+=3;
	  }
      }
    clean_vector(dataI_in);
    vector<int>head_number;
    for(int ni=0;ni<NodeGroups.size();ni++)
      head_number.push_back(ni);
    for(int i=0; i < list_pair_1.size(); ++i)
      {
	int j=list_pair_1[i];
	while(head_number[j] != j)
	  j=head_number[j];
	int k=list_pair_2[i];
	while(head_number[k] != k)
	  k=head_number[k];
	if(j != k) 
	  head_number[j]=k;
      }
    clean_vector(list_pair_1);
    clean_vector(list_pair_2);
    for (int j=0;j < head_number.size(); ++j)
      while(head_number[j] != head_number[head_number[j]])
	head_number[j]=head_number[head_number[j]];
    set<int>heads;
    for (int j=0;j < head_number.size(); ++j)
      if(head_number[j]==j)
	heads.insert(j);
    vector<map<pair<int,int>,int>>tmpSGroups(heads.size());
    vector<int>TotalPoints(heads.size(),0);
    int j=0;
    for(auto NG : NodeGroups)
      {
	auto it=heads.find(head_number[j++]);
	int k=distance(heads.begin(),it);
	tmpSGroups[k].insert(NG);
	TotalPoints[k]+=NG.second;
      }
    clean_vector(head_number);
    heads.clear();
    multimap<int,int>Tpoints;
    int TPcount=0;
    for(auto TP : TotalPoints)
      Tpoints.insert({-TP,TPcount++});
    deque<map<pair<int,int>,int>>SGroups;
    for(auto TP : Tpoints)
      SGroups.push_back(tmpSGroups[TP.second]);
    tmpSGroups.clear();
    deque<set<int>> NODESinSG;
    int countc=0;
    for(auto SG : SGroups)
      {
	NODESinSG.push_back(set<int>());
	int countb=0;
	for(auto sg : SG)
	  NODESinSG.back().insert(sg.first.first);
	for(auto FR : NODESinSG.back())
	  {
	    FHT << " INDEP " << mem.steps << " " << level << " " << countb++ << " " << countc << FR << "\n";
	  }
	countc++;
      }
    while (!SGroups.empty())
      {
	set<int>FRused;
	S2NGroups.push_back(vector<map<pair<int,int>,int>>());
	S2Nodes.push_back(vector<set<int>>());
	FreeNodes.push_back(vector<int>());
	auto itSG=SGroups.begin();
	auto itFR=NODESinSG.begin();
	while(itFR!=NODESinSG.end())
	  {
	    set<int>FRtmp;
	    for(auto FR : *itFR)
	      {
		if(FRused.find(FR) != FRused.end())
		  {
		    FRtmp.clear();
		    break;
		  }
		else
		  FRtmp.insert(FR);
	      }
	    if(!FRtmp.empty())
	      {
		S2NGroups.back().push_back(*itSG);
		S2Nodes.back().push_back(*itFR);
		SGroups.erase(itSG);
		NODESinSG.erase(itFR);
		FRused.insert(FRtmp.begin(),FRtmp.end());
	      }
	    else
	      {
		itSG++;
		itFR++;
	      }
	  }
      }
    for(int W=1;W < S2NGroups.size();W++)
      {
	WorldRanks.push_back(vector<int>());
	LocalGroups.push_back(vector<int>());
	FreeNodes.push_back(vector<int>());
	IAmIn.push_back(bool());
	vector<bool>useit(FractalNodes,false);
	int count=-1;
	for(auto SN : S2Nodes[W])
	  {
	    count++;
	    for(int N : SN)
	      useit[N]=true;
	    IAmIn.back()= SN.find(FractalRank) != SN.end();
	    if(!IAmIn.back())
	      continue;
	    WorldRanks.back().resize(SN.size());
	    copy(SN.begin(),SN.end(),WorldRanks.back().begin());
	    for(auto NG : S2NGroups[W][count])
	      if(NG.first.first == FractalRank)
		LocalGroups.back().push_back(NG.first.second);
	    break;
	  }
	for(int FR=0;FR<FractalNodes;FR++)
	  if(!useit[FR])
	    FreeNodes.back().push_back(FR);
      }
  }
}
