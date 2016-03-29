#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void add_buffer_values(Fractal_Memory& mem,int level,
			 vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    int FractalRank=mem.p_mess->FractalRank;
    vector <int>BOX=mem.BoxesLev[FractalRank][level];
    int HypreNodes=mem.p_mess->HypreNodes;
    vector <int> counts_in(HypreNodes);
    vector <int> counts_out(HypreNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(HypreNodes);
    vector < vector <double> > dataR_out(HypreNodes);
    vector <int>pos(3);
    // for(int TW=0;TW<mem.Touchy.size();TW++)
    for(int FR : mem.Touchy)
      {
	// int FR=mem.Touchy[TW];
	int HR=mem.p_mess->IHranks[FR];
	assert(HR >= 0);
	vector <int>TBBOX=mem.BBoxesLev[FR][level];
	// for(int B=0;B<SBoxes.size();B++)
	int B=0;
	for(vector <int>& SB : SBoxes)
	  {
	    if(overlap_boxes(TBBOX,SB))
	      {
		for(Point* &p : SPoints[B])
		  {
		    p->get_pos_point(pos);
		    if(vector_in_box(pos,TBBOX))
		      {
			dataI_out[HR].push_back(pos[0]);
			dataI_out[HR].push_back(pos[1]);
			dataI_out[HR].push_back(pos[2]);
			dataR_out[HR].push_back(p->get_potential_point());
			counts_out[HR]++;
		      }
		  }
	      }
	    B++;
	  }
      }
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=3;
    int doubles=1;
    mem.p_mess->Send_Data_Some_How(9,mem.p_mess->HypreWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    std::set<Point*>edgeP;
    for(Group* pg : mem.all_groups[level])
      {
	for(Point* &p : pg->list_points)
	  {
	    if(!p->get_inside() && !vector_in_box(pos,BOX))
	      edgeP.insert(p);
	  }
      }
    auto iend=edgeP.end();
    int c1=0;
    int c3=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	for(int c=0;c<counts_in[HR];c++)
	  {
	    Point* p;
	    p->set_pos_point(dataI_in[c3],dataI_in[c3+1],dataI_in[c3+2]);
	    auto it=edgeP.find(p);
	    if(it != iend)
	      (*it)->set_potential_point(dataR_in[c1]);
	    c1++;
	    c3+=3;
	  }
      }
  }
}
