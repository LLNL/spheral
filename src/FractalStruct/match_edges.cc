#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void match_edges(Fractal_Memory& mem,int level)
  {
    ofstream& FHT=mem.p_file->DUMPS;
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <int> counts_in(FractalNodes);
    vector <int> counts_out(FractalNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(FractalNodes);
    vector < vector <double> > dataR_out(FractalNodes);
    map <array<int,3>,Point*,point_comp2> lows;
    int nna=0;
    int nnb=0;
    int nnc=0;
    int nnd=0;
    int nne=0;
    int sumta=0;
    int sumtb=0;
    for(auto pg : mem.all_groups[level])
      {
	// FHT << " MEA " << nna << endl;
	for(auto &p : pg->list_points)
	  {
	    // FHT << " MEB " << nnb << endl;
	    if(p->get_edge_point() || p->get_buffer_point())
	      {
		// FHT << " MEC " << nnc << endl;
		array<int,3>pos=p->get_pos_point_a();
		if(p->get_it_is_high())
		  {
		    // FHT << " MED " << nnd << endl;
		    for(int FR : mem.TouchWhichBoxes)
		      {
			if(!vector_in_box(pos,mem.BBoxesLev[FR][level]))
			  continue;
			dataI_out[FR].push_back(pos[0]);
			dataI_out[FR].push_back(pos[1]);
			dataI_out[FR].push_back(pos[2]);
			counts_out[FR]++;
		      }
		    nnd++;
		  }
		else
		  {
		    FHT << " MEE " << nne << endl;
		    lows[pos]=p;
		    nne++;
		  }
		nnc++;
	      }
	    nnb++;
	  }
	nna++;
      }
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=3;
    int doubles=0;
    FHT << " MEF " << nna << " " << nnb << " " << nnc << " " << nnd << " " << nne << " " << endl;
    for(int FR=0;FR<FractalNodes;FR++)
      FHT << "OUTME " << level << " " << FR << " " << dataI_out[FR].size() << endl;
    cerr << " MATCHA " << FractalRank << endl;
    mem.p_mess->Full_Stop_Do_Not_Argue();
    mem.p_mess->Send_Data_Some_How(39,mem.p_mess->FractalWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    cerr << " MATCHB " << FractalRank << endl;
    FHT << " MEG " << endl;
    mem.p_mess->Full_Stop_Do_Not_Argue();
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    int c3=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    array<int,3>ar3={dataI_in[c3],dataI_in[c3+1],dataI_in[c3+2]};
	    c3++;
	    auto f=lows.find(ar3);
	    if(f != lows.end())
	      {
		f->second->set_it_is_high(true);
		FHT << " MATCH " << level << " " << FR << " " << ar3[0] << " " << ar3[1] << " " << ar3[2] << " " << "\n";
	      }
	  }
      }
  }
}
