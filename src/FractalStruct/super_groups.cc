#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void super_groups(Fractal_Memory& mem,vector <Group*>& groups,const int level)
  {
    const int FractalNodes=mem.p_mess->FractalNodes;
    const int FractalRank=mem.p_mess->FractalRank;
    const int ROOTNODE=mem.p_mess->ROOTNODE;
    int spacing=Misc::pow(2,mem.level_max-level);
    vector<int>LowTouch;
    vector<vector<int>>LowBox;
    vector<int>HighTouch;
    vector<vector<int>>HighBox;
    vector <int>MyBox;
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
	else
	  {
	    MyBox=mem.BoxesLev[FractalRank][level];
	    for(int ni : {1,3,5})
	      MyBox[ni]+=spacing;
	  }
      }
    map<pair<int,int>,int>NodeGroups;
    { // Braces are here deliberately.
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
      int group_number=-1;
      for(auto pgroup : groups)
	{
	  group_number++;
	  if(!pgroup->get_buffer_group())
	    continue;
	  dataI_out[ROOTNODE].push_back(group_number);
	  dataI_out[ROOTNODE].push_back(pgroup->list_points.size());
	  counts_out[ROOTNODE]++;
	}
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
    std::map<array<int,3>,int,point_comp2>SposG;
    std::set<array<int,3>,point_comp2> GroupHug;
    { // Braces are here deliberately.
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
	    SposG.insert(make_pair(pos,group_number));// ignore overlap for now.
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
	    auto it=SposG.find({{dataI_in[c4],dataI_in[c4+1],dataI_in[c4+2]}});
	    if(it != SposG.end())
	      GroupHug.insert({{it->second,FR,dataI_in[c4+3]}});
	    c4+=4;
	  }
    }
    vector<int>head_number;
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
      for(auto GH : GroupHug)
	{
	  dataI_out[ROOTNODE].push_back(GH[0]);
	  dataI_out[ROOTNODE].push_back(GH[1]);
	  dataI_out[ROOTNODE].push_back(GH[2]);
	  counts_out[ROOTNODE]++;
	}
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
    clean_vector(list_pair_1);
    clean_vector(list_pair_2);
    for (int j=0;j < NodeGroups.size(); ++j)
      {
	while(head_number[j] != head_number[head_number[j]])
	  head_number[j]=head_number[head_number[j]];
      }
    set<int>heads;
    for (int j=0;j < NodeGroups.size(); ++j)
      if(head_number[j]==j)
	heads.insert(j);
    vector<vector<int>>SGroups(heads.size());
    vector<int>TotalPoints(heads.size(),0);
    for (int j=0;j < NodeGroups.size(); ++j)
      {
	auto it=heads.find(head_number[j]);
	int k=distance(heads.begin(),it);
	SGroups[k].push_back(j);
	// auto itj=next(NodeGroups.begin(),j);
	TotalPoints[k]+=next(NodeGroups.begin())->second;
      }
    clean_vector(head_number);
    heads.clear();
    ofstream& FHT=mem.p_file->DUMPS;
    for(int ni=0;ni<heads.size();ni++)
      {
	FHT << " SUPER " << mem.steps << " " << mem.level << " " << ni << " " << TotalPoints[ni] << " " << SGroups[ni].size() << "\n";
      }
  }
}
