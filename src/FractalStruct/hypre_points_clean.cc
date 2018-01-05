#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_clean(Fractal_Memory& mem,int level,vector< vector<Point*> >& hypre_points)
  {
    ofstream& FHT=mem.p_file->DUMPS;
    int FractalRank=mem.p_mess->FractalRank;      
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <int>Touchies;
    for(int FR : mem.TouchWhichBoxes)
      if(mem.p_mess->counts_on_nodes[2*level+1][FR])
	{
	  Touchies.push_back(FR);
	  FHT << " Touchies " << FR << endl;
	}
    map<array<int,3>,set<int>,point_comp2> edgies;
    map<array<int,3>,Point*,point_comp2> mepoints;
    vector <int>BOX=mem.BoxesLev[FractalRank][level];
    {
      vector <int> counts_in(FractalNodes);
      vector <int> dataI_in;
      {
	vector <int> counts_out(FractalNodes,0);
	vector <double> dataR_in;
	vector <vector<int>> dataI_out(FractalNodes);
	vector <vector<double>> dataR_out(FractalNodes);
	for(auto &ph : hypre_points)
	  for(auto &p : ph)
	    {
	      if(!on_edge(p,BOX))
		continue;
	      array<int,3>epos(p->get_pos_point_a());
	      mepoints.insert(make_pair(epos,p));
	      set<int>dirs;
	      for(int ni=0;ni<6;ni++)
		{
		  Point* p1=p->get_point_ud_0(ni,-13);
		  if(vector_in_box(p1,BOX))
		    dirs.insert(ni);
		  else
		    {
		      array<int,3>pos1(p1->get_pos_point_a());
		      int count(0);
		      for(int FR : Touchies)
			{
			  if(!on_edge(pos1,mem.BoxesLev[FR][level]))
			    continue;
			  dataI_out[FR].push_back(pos1[0]);
			  dataI_out[FR].push_back(pos1[1]);
			  dataI_out[FR].push_back(pos1[2]);
			  int nia=(ni%2) == 0 ? ni+1:ni-1;
			  dataI_out[FR].push_back(nia);
			  counts_out[FR]++;
			  assert((++count) < 2);
			}
		    }
		}
	      array<int,3>pos(p->get_pos_point_a());
	      auto ret=edgies.insert(make_pair(pos,dirs));
	      assert(ret.second);
	    }
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=4;
	int doubles=0;
	mem.p_mess->Send_Data_Some_How(12,mem.p_mess->FractalWorld,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
      }
      int c4(0);
      for(int FR=0;FR<FractalNodes;FR++)
	{
	  for(int c=0;c<counts_in[FR];c++)
	    {
	      auto it=edgies.find(array<int,3>{{dataI_in[c4],dataI_in[c4+1],dataI_in[c4+2]}});
	      if(it != edgies.end())
		it->second.insert(dataI_in[c4+3]);
	      c4+=4;
	    }
	}
    }
    vector<int>totals(10,0);
    int maxy=0;
    for(auto &p : edgies)
      {
	int what=p.second.size();
	totals[what]++;
	maxy=max(maxy,what);
	if(what != 6)
	  {
	    array<int,3>pos=p.first;
	    auto itme=mepoints.find(pos);
	    for(int ni=0;ni<6;ni++)
	      {
		if(p.second.find(ni) == p.second.end())
		  {
		    Point* tr=itme->second->get_point_ud_0(ni,-15);
		    tr->set_trouble(true);
		    // array<int,3>posa=tr->get_pos_point_a();
		    // FHT << " EPA " << BOX[0] << " " << posa[0] << " " << BOX[1] << "   "
		    // 	<< BOX[2] << " " << posa[1] << " " << BOX[3] << "   "
		    // 	<< BOX[4] << " " << posa[2] << " " << BOX[5] << " " << ni << "\n";
		  }
	      }
	    // bool down=pos[0] == BOX[0] || pos[1] == BOX[2] || pos[2] == BOX[4];
	    // bool up=pos[0] == BOX[1] || pos[1] == BOX[3] || pos[2] == BOX[5];
	    // string d = down ? " D " : " d ";
	    // string u = up ? " U " : " u ";
	    // FHT << " EPB " << BOX[0] << " " << pos[0] << " " << BOX[1] << "   "
	    // 	<< BOX[2] << " " << pos[1] << " " << BOX[3] << "   "
	    // 	<< BOX[4] << " " << pos[2] << " " << BOX[5] << d << u << what << " ""\n";
	  }
      }
    for(int t=0;t<=maxy;t++)
      FHT << " EDGES " << level << " " << t << " " << totals[t] << " " << edgies.size() << "\n";
  }
}
