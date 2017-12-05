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
    WorldRanks.clear();
    LocalGroups.clear();
    FreeNodes.clear();
    clean_vector(IAmIn);
    clean_vector(mem.Touchy);
    WorldRanks.resize(1);
    LocalGroups.resize(1);
    FreeNodes.resize(1);
    IAmIn.resize(1,false);
    const int FractalNodes=mem.p_mess->FractalNodes;
    const int FractalRank=mem.p_mess->FractalRank;
    int spacing=Misc::pow(2,mem.level_max-level);
    vector<int>LowTouch;
    vector<vector<int>>LowBox;
    vector<int>HighTouch;
    vector<vector<int>>HighBox;
    vector <int>MyBox(mem.BoxesLev[FractalRank][level]);
    for(int ni : {1,3,5})
      MyBox[ni]+=spacing;
    
    int counta=0;
    int countb=0;
    for(auto pg : groups)
      if(pg->get_buffer_group())
	countb=2;
      else
	counta=1;
    vector<int>counts(FractalNodes);
    mem.p_mess->How_Many_On_Nodes(counta+countb,counts);
    vector<bool>occupieda(FractalNodes);
    vector<bool>occupiedb(FractalNodes);
    for(int FR=0;FR<FractalNodes;FR++)
      {
	occupieda[FR]=counts[FR] % 2 == 1;
	occupiedb[FR]=(counts[FR]/2) % 2 == 1;
      }
    for(auto FR : mem.TouchWhichBoxes)
      {
	if(counts[FR] == 0 || counts[FractalRank] == 0)
	  continue;
	vector <int>TBBox=mem.BBoxesLev[FR][level];
	for(auto pg : groups)
	  {
	    if(group_in_box(pg,TBBox))
	      {
		mem.Touchy.push_back(FR);
		break;
	      }
	  }
      }
    for(int FR : mem.Touchy)
      {
	if(FR < mem.p_mess->FractalRank)
	  {
	    LowTouch.push_back(FR);
	    LowBox.push_back(mem.BoxesLev[FR][level]);
	    for(int ni : {1,3,5})
	      LowBox.back()[ni]+=spacing;
	  }
	else if(FR > mem.p_mess->FractalRank)
	  {
	    HighTouch.push_back(FR);
	    HighBox.push_back(mem.BoxesLev[FR][level]);
	    for(int ni : {1,3,5})
	      HighBox.back()[ni]+=spacing;
	  }
      }
    vector<vector<map<pair<int,int>,int>>>S2NGroups(1);
    vector<vector<set<int>>>S2Nodes(1);
    map<pair<int,int>,int>NodeGroups; // FR, Group#, Points in Group
    { // Braces are here deliberately.
      int group_number=-1;
      for(auto pgroup : groups)
	{
	  group_number++;
	  if(pgroup->get_buffer_group())
	    continue;
	  WorldRanks.back().resize(1);
	  WorldRanks.back().back()=FractalRank;
	  LocalGroups.back().push_back(group_number);
	  IAmIn.back()=true;
	  S2NGroups.back().resize(1);
	  S2Nodes.back().resize(1);
	  S2NGroups.back().back().insert(make_pair(make_pair(FractalRank,group_number),pgroup->list_points.size()));
	  S2Nodes.back().back().insert(FractalRank);
	}
      for(int FR=0;FR<FractalNodes;FR++)
	if(!occupieda[FR])
	  FreeNodes.back().push_back(FR);
      vector <int> counts_in(FractalNodes);
      vector <int> counts_out(FractalNodes,0);
      vector <int> dataI_in;
      vector <double> dataR_in;
      vector < vector <int> > dataI_out(FractalNodes);
      vector < vector <double> > dataR_out(FractalNodes);
      int how_manyI=-1;
      int how_manyR=-1;
      int integers=2;
      int doubles=0;
      group_number=-1;
      vector<int>DI;
      for(auto pgroup : groups)
	{
	  group_number++;
	  if(!pgroup->get_buffer_group())
	    continue;
	  DI.push_back(group_number);
	  DI.push_back(pgroup->list_points.size());
	}
      const int co=DI.size()/2;
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
      int c2=0;
      for(int FR=0;FR<FractalNodes;FR++)
	for(int c=0;c<counts_in[FR];c++)
	  {
	    NodeGroups[make_pair(FR,dataI_in[c2])]=dataI_in[c2+1];
	    c2+=2;
	  }
    }
    std::set<array<int,3>,point_comp2> GroupHug;
    { // Braces are here deliberately.
      std::map<array<int,3>,int,point_comp2>MposG;
      vector <int> counts_in(FractalNodes);
      vector <int> counts_out(FractalNodes,0);
      vector <int> dataI_in;
      vector <double> dataR_in;
      vector < vector <int> > dataI_out(FractalNodes);
      vector < vector <double> > dataR_out(FractalNodes);
      int how_manyI=-1;
      int how_manyR=-1;
      int integers=4;
      int doubles=0;
      int group_number=-1;
      for(auto pgroup : groups)
	{
	  group_number++;
	  if(!pgroup->get_buffer_group())
	    continue;
	  for(Point* p : pgroup->list_points)
	    {
	      if(!on_edge(p,MyBox))
		continue;
	      array<int,3>pos=p->get_pos_point_a();
	      MposG.insert(make_pair(pos,group_number));// ignore overlap for now.
	      int LB=-1;
	      for(int FR : LowTouch)
		{
		  LB++;
		  if(!on_edge(pos,LowBox[LB]))
		    continue;
		  dataI_out[FR].push_back(pos[0]);
		  dataI_out[FR].push_back(pos[1]);
		  dataI_out[FR].push_back(pos[2]);
		  dataI_out[FR].push_back(group_number);
		  counts_out[FR]++;
		}
	    }	
	}
      mem.p_mess->Send_Data_Some_How(29,mem.p_mess->FractalWorld,counts_out,counts_in,integers,doubles,
				     dataI_out,dataI_in,how_manyI,
				     dataR_out,dataR_in,how_manyR);
      clean_vector(counts_out);
      dataI_out.clear();
      dataR_out.clear();
      int c4=0;
      for(int FR=0;FR<FractalNodes;FR++)
	for(int c=0;c<counts_in[FR];c++)
	  {
	    auto it=MposG.find({{dataI_in[c4],dataI_in[c4+1],dataI_in[c4+2]}});
	    if(it != MposG.end())
	      GroupHug.insert({{it->second,FR,dataI_in[c4+3]}});
	    c4+=4;
	  }
    }
    vector<int>head_number;
    {
      vector<int>list_pair_1;
      vector<int>list_pair_2;
      {
	vector <int> counts_in(FractalNodes);
	vector <int> counts_out(FractalNodes,0);
	vector <int> dataI_in;
	vector <double> dataR_in;
	vector < vector <int> > dataI_out(FractalNodes);
	vector < vector <double> > dataR_out(FractalNodes);
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=3;
	int doubles=0;
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
	for(int ni=0;ni<NodeGroups.size();ni++)
	  head_number.push_back(ni);
	int c3=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		auto it=NodeGroups.find(make_pair(FR,dataI_in[c3]));
		list_pair_1.push_back(distance(NodeGroups.begin(),it));
		it=NodeGroups.find(make_pair(dataI_in[c3+1],dataI_in[c3+2]));
		list_pair_2.push_back(distance(NodeGroups.begin(),it));
	      }
	    c3+=3;
	  }
      }
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
      for (int j=0;j < NodeGroups.size(); ++j)
	while(head_number[j] != head_number[head_number[j]])
	  head_number[j]=head_number[head_number[j]];
    }
    set<int>heads;
    for (int j=0;j < NodeGroups.size(); ++j)
      if(head_number[j]==j)
	heads.insert(j);
    vector<map<pair<int,int>,int>>tmpSGroups(heads.size());
    vector<int>TotalPoints(heads.size(),0);
    int j=0;
    for(auto NG : NodeGroups)
      {
	auto it=heads.find(head_number[j]);
	int k=distance(heads.begin(),it);
	tmpSGroups[k].insert(NG);
	TotalPoints[k]+=NG.second;
	j++;
      }
    clean_vector(head_number);
    heads.clear();
    multimap<int,int>Tpoints;
    int TPcount=0;
    for(auto TP : TotalPoints)
      Tpoints.insert(make_pair(-TP,TPcount++));
    deque<map<pair<int,int>,int>>SGroups;
    for(auto TP : Tpoints)
      SGroups.push_back(tmpSGroups[TP.second]);
    tmpSGroups.clear();
    deque<set<int>> NODESinSG;
    for(auto SG : SGroups)
      {
	set<int>FRtmp;
	for(auto sg : SG)
	  FRtmp.insert(std::get<0>(sg.first));
	NODESinSG.resize(NODESinSG.size()+1);
	for(auto FR : FRtmp)
	  NODESinSG.back().insert(FR);
      }
    while (!SGroups.empty())
      {
	set<int>FRused;
	S2NGroups.resize(S2NGroups.size()+1);
	S2Nodes.resize(S2Nodes.size()+1);
	FreeNodes.resize(FreeNodes.size()+1);
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
    int Worlds=S2NGroups.size();
    for(int W=1;W<Worlds;W++)
      {
	WorldRanks.resize(WorldRanks.size()+1);
	LocalGroups.resize(LocalGroups.size()+1);
	FreeNodes.resize(FreeNodes.size()+1);
	IAmIn.resize(IAmIn.size()+1);
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
	      {
		auto p=NG.first;
		if(get<0>(p) == FractalRank)
		  LocalGroups.back().push_back(get<1>(p));
	      }
	    break;
	  }
	for(int FR=0;FR<FractalNodes;FR++)
	  if(!useit[FR])
	    FreeNodes.back().push_back(FR);
      }
  }
}
