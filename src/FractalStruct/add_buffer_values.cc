#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void add_buffer_values(Fractal_Memory& mem,int level,
			 vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    int HypreNodes=mem.p_mess->HypreNodes;
    int spacing=Misc::pow(2,mem.p_fractal->get_level_max()-level);
    // for(int B=0;B<SBoxes.size();B++)
    //   Misc::times(SBoxes[B],spacing);
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
	if(HR < 0 )
	  continue;
	vector <int>TBBOX=mem.BBoxesLev[FR][level];
	// for(int B=0;B<SBoxes.size();B++)
	int B=0;
	for(vector <int>& SB : SBoxes)
	  {
	    if(overlap_boxes(TBBOX,SB))
	      {
		dataI_out[HR].push_back(B);
		dataI_out[HR].insert(dataI_out[HR].end(),SB.begin(),SB.end());
		counts_out[HR]++;
	      }
	    B++;
	  }
      }
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=7;
    int doubles=0;
    mem.p_mess->Send_Data_Some_How(9,mem.p_mess->HypreWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();

    dataI_out.resize(HypreNodes);
    dataR_out.resize(HypreNodes);
    counts_out.assign(HypreNodes,0);
    int counterI=0;
    auto pcounter=dataI_in.begin();
    for(int HR=0;HR<HypreNodes;HR++)
      {
	int FR=mem.p_mess->Hranks[HR];
	vector <int>TBBOX;
	vector <int>Box_coords(6);
	vector <int>Blist;
	for(int c=0;c<counts_in[HR];c++)
	  {
	    if(c == 0)
	      {
		TBBOX=mem.BBoxesLev[FR][level];
		Blist.clear();
		// for(int B=0;B<SBoxes.size();B++)
		int B=0;
		for(auto& SB : SBoxes)
		  if(overlap_boxes(TBBOX,SB))
		    Blist.push_back(B);
		B++;
	      }
	    int Box_number=*pcounter;
	    std::copy(pcounter+1,pcounter+7,Box_coords.begin());
	    for(int B=0;B<Blist.size();B++)
	      {
		if(!overlap_boxes(Box_coords,SBoxes[Blist[B]]))
		  continue;
		for(vector <Point*>::const_iterator p_itr=SPoints[Blist[B]].begin();p_itr!=SPoints[Blist[B]].end();p_itr++)
		  {
		    Point* p=*p_itr;
		    p->get_pos_point(pos);
		    if(!vector_in_box(pos,TBBOX))
		      continue;
		    for(int ni=0;ni<6;ni++)
		      {
			Point* p1=p->get_point_ud(ni);
			p1->get_pos_point(pos);
			if(!vector_in_box(pos,Box_coords))
			  continue;
			dataI_out[HR].push_back(Box_number);
			dataI_out[HR].push_back(6*Misc::coordinate(pos,Box_coords,spacing) + (ni % 2 == 0 ? ni+1:ni-1));
			dataR_out[HR].push_back(p1->get_potential_point());
			counts_out[HR]++;
			break;
		      }
		  }
	      }
	    std::advance(pcounter,7);
	  }
      }
    dataI_in.clear();
    dataR_in.clear();
    counts_in.assign(HypreNodes,0);
    integers=2;
    doubles=1;
    how_manyI=-1;
    how_manyR=-1;
    mem.p_mess->Send_Data_Some_How(10,mem.p_mess->HypreWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    counterI=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	for(int c=0;c<counts_in[HR];c++)
	  {
	    int B=dataI_in[counterI];
	    int number=dataI_in[counterI+1]/6;
	    int dir=dataI_in[counterI+1] % 6;
	    SPoints[B][number]->get_point_ud(dir)->set_potential_point(dataR_in[counterI/2]);
	    counterI+=2;
	  }
      }
    for(int B=0;B<SBoxes.size();B++)
      Misc::divide(SBoxes[B],spacing);
  }
}
