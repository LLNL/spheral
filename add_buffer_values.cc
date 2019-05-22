#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void add_buffer_values(Fractal_Memory& mem,int level)
  {
    static int _COUNTER=0;
    ofstream& FF = mem.p_file->DUMPS;
    const int FractalRank=mem.p_mess->FractalRank;
    // const int FractalNodes=mem.p_mess->FractalNodes;
    // const int HypreRank=mem.p_mess->HypreRank;
    const int HypreNodes=mem.p_mess->HypreNodes;
    vector <int> counts_in(HypreNodes);
    vector <int> counts_out(HypreNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(HypreNodes);
    vector < vector <double> > dataR_out(HypreNodes);
    vector<int>Box=mem.BoxesLev[FractalRank][level];
    FF << " ALLA " << mem.steps << " " << Box[0] << " " << Box[1] << " " << Box[2] << " " << Box[3] << " " << Box[4] << " " << Box[5] << "\n";
    vector<int>BBox=mem.BBoxesLev[FractalRank][level];
    FF << " ALLA " << mem.steps << " " << BBox[0] << " " << BBox[1] << " " << BBox[2] << " " << BBox[3] << " " << BBox[4] << " " << BBox[5] << "\n";
    int WIDTH=mem.periodic ? mem.grid_length*Misc::pow(2,mem.level_max) : 0;
    std::map<std::array<int,3>,Point*,point_comp2> edgeP;
    vector<Point*>bvec;
    int efound=0;
    int ifound=0;
    for(Group* &pg : mem.all_groups[level])
      {
	if(!pg->get_buffer_group())
	  continue;
	for(Point* &p : pg->list_points)
	  {
	    double pott=p->get_potential_point();
	    auto ar3=p->get_pos_point_a();
	    if(p->get_edge_point())
	      {
		auto succ=edgeP.insert(make_pair(ar3,p));
		// FF << " FOUNDEDGE " << mem.steps << " ";
		// FF << efound++ << " " << ifound << " " << ar3[0] << " " << ar3[1] << " " << ar3[2]  << " " << pott << " " << p <<"\n";
		if(!succ.second)
		  {
		    ar3=succ.first->first;
		    FF << " FOUNDDUP " << mem.steps << " ";
		    FF << efound << " " << ifound << " " << ar3[0] << " " << ar3[1] << " " << ar3[2]  << " " << pott << " " << succ.first->second << "\n";

		  }
	      }
	    else if(on_edge(ar3,Box))
	      {
		bvec.push_back(p);
		// FF << " FOUNDIN " << mem.steps << " ";
		// FF << efound << " " << ifound++ << " " << ar3[0] << " " << ar3[1] << " " << ar3[2]  <<" " <<  pott << "\n";
	      }
	  }
      }

    for(int HR=0;HR<HypreNodes;HR++)
      {
	int FR=mem.p_mess->Hranks[HR];
	if(FR == FractalRank)
	  continue;
	vector <int>FRBBOX=mem.BBoxesLev[FR][level];
	for(auto p : bvec)
	  {
	    vector <int>pos=p->get_pos_point();
	    if(vector_in_box(pos,FRBBOX,WIDTH,mem.periodic))
	      {
		dataI_out[HR].push_back(pos[0]);
		dataI_out[HR].push_back(pos[1]);
		dataI_out[HR].push_back(pos[2]);
		dataR_out[HR].push_back(p->get_potential_point());
		counts_out[HR]++;
	      }
	  }
      }
    clean_vector(bvec);
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=3;
    int doubles=1;
    mem.p_mess->Send_Data_Some_How(9,mem.p_mess->HypreWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    clean_vector(counts_out);
    dataI_out.clear();
    dataR_out.clear();

    int c1=0;
    int c3=0;
    int found=0;
    int notfound=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	for(int c=0;c<counts_in[HR];c++)
	  {
	    array <int,3>posin={dataI_in[c3],dataI_in[c3+1],dataI_in[c3+2]};
	    auto eb=edgeP.find(posin);
	    if(eb != edgeP.end())
	      {
		// double pott=eb->second->get_potential_point();
		eb->second->set_potential_point(dataR_in[c1]);
		// FF << " FOUNDIT " << mem.steps << " " << mem.p_mess->Hranks[HR] << " ";
		// FF << found+1 << " " << notfound << " " << posin[0] << " " << posin[1] << " " << posin[2] << " " << dataR_in[c1] << " " << pott << "\n";
		found++;
	      }
	    else
	      {
		// FF << " FOUNDNO " << mem.steps << " " << mem.p_mess->Hranks[HR] << " ";
		// FF << found << " " << notfound+1 << " " << posin[0] << " " << posin[1] << " " << posin[2] << " " << dataR_in[c1] << "\n";
		notfound++;
	      }
	    c1++;
	    c3+=3;
	  }
      }
    _COUNTER++;
  }
}
