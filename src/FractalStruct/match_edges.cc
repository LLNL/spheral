#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void match_edges(Fractal_Memory& mem,int level)
  {
    // ofstream& FHT=mem.p_file->DUMPS;
    // int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <int> counts_in(FractalNodes);
    vector <int> counts_out(FractalNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(FractalNodes);
    vector < vector <double> > dataR_out(FractalNodes);
    multimap <array<int,3>,Point*,point_comp2> lows;
    array<int,3>ar3;
    int npg=0;
    int np=0;
    for(auto &pg : mem.all_groups[level])
      {
	for(auto &p : pg->list_points)
	  {
	    if(p->get_buffer_point() && !p->get_it_is_high())
	      {
		array<int,3>ar3=p->get_pos_point_a();
		// FHT << " AR3 " << level << " " << npg << " " << np << " " << ar3[0] << " " << ar3[1] << " " << ar3[2] << " " << p << endl;
		lows.insert(pair<array<int,3>,Point*>(ar3,p));
	      }
	    if(p->get_edge_point() && p->get_it_is_high() && !p->get_it_is_really_high())
	      {
		for(int FR : mem.TouchWhichBoxes)
		  {
		    if(!vector_in_box(ar3,mem.BBoxesLev[FR][level]))
		      continue;
		    dataI_out[FR].push_back(ar3[0]);
		    dataI_out[FR].push_back(ar3[1]);
		    dataI_out[FR].push_back(ar3[2]);
		    counts_out[FR]++;
		  }
	      }
	    np++;
	  }
	npg++;
      }
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=3;
    int doubles=0;
    mem.p_mess->Send_Data_Some_How(39,mem.p_mess->FractalWorld,
				   counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    int c3=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    array<int,3>ar3{{dataI_in[c3],dataI_in[c3+1],dataI_in[c3+2]}};
	    c3++;
	    auto range=lows.equal_range(ar3);
	    for(auto found=range.first;found!=range.second;found++)
	      found->second->set_it_is_high(true);
	  }
      }
  }
}
