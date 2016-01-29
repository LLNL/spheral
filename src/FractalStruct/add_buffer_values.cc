#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void add_buffer_values(Fractal_Memory& mem,int level,
			 vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    int FractalNodes=mem.FractalNodes;
    int spacing=Misc::pow(2,mem.p_fractal->get_level_max()-level);
    vector <int> counts_in(FractalNodes);
    vector <int> counts_out(FractalNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(FractalNodes);
    vector < vector <double> > dataR_out(FractalNodes);
    vector <int>pos(3);
    for(int TW=0;TW<mem.Touchy.size();TW++)
      {
	int FR=mem.Touchy[TW];
	vector <int>TBBOX=mem.BBoxesLev[FR][level];
	for(int B=0;B<SBoxes.size();B++)
	  {
	    if(overlap_boxes(TBBOX,SBoxes[B]))
	      {
		dataI_out[FR].push_back(B);
		dataI_out[FR].push_back(SBoxes[B][0]);
		dataI_out[FR].push_back(SBoxes[B][1]);
		dataI_out[FR].push_back(SBoxes[B][2]);
		dataI_out[FR].push_back(SBoxes[B][3]);
		dataI_out[FR].push_back(SBoxes[B][4]);
		dataI_out[FR].push_back(SBoxes[B][5]);
		counts_out[FR]++;
	      }
	  }
      }
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=7;
    int doubles=0;
    mem.p_mess->Send_Data_Some_How(9,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();

    dataI_out.resize(FractalNodes);
    dataR_out.resize(FractalNodes);
    counts_out.assign(FractalNodes,0);
    int counterI=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	vector <int>TBBOX;
	vector <int>Box_coords(6);
	vector <int>Blist;
	for(int c=0;c<counts_in[FR];c++)
	  {
	    if(c == 0)
	      {
		TBBOX=mem.BBoxesLev[FR][level];
		Blist.clear();
		for(int B=0;B<SBoxes.size();B++)
		  {
		    if(overlap_boxes(TBBOX,SBoxes[B]))
		      Blist.push_back(B);
		  }
	      }
	    int Box_number=dataI_in[counterI];
	    Box_coords[0]=dataI_in[counterI+1];
	    Box_coords[1]=dataI_in[counterI+2];
	    Box_coords[2]=dataI_in[counterI+3];
	    Box_coords[3]=dataI_in[counterI+4];
	    Box_coords[4]=dataI_in[counterI+5];
	    Box_coords[5]=dataI_in[counterI+6];
	    for(int B=0;B<Blist.size();B++)
	      {
		if(!overlap_boxes(Box_coords,SBoxes[Blist[B]]))
		  continue;
		for(vector <Point*>::const_iterator p_itr=SPoints[Blist[B]].begin();p_itr<SPoints[Blist[B]].end();p_itr++)
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
			dataI_out[FR].push_back(Box_number);
			dataI_out[FR].push_back(6*Misc::coordinate(pos,Box_coords,spacing) + (ni % 2 == 0 ? ni+1:ni-1));
			dataR_out[FR].push_back(p1->get_potential_point());
			counts_out[FR]++;
			break;
		      }
		  }
	      }
	    counterI+=7;
	  }
      }
    dataI_in.clear();
    dataR_in.clear();
    counts_in.assign(FractalNodes,0);
    integers=2;
    doubles=1;
    how_manyI=-1;
    how_manyR=-1;
    mem.p_mess->Send_Data_Some_How(10,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    counterI=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    int B=dataI_in[counterI];
	    int number=dataI_in[counterI+1]/6;
	    int dir=dataI_in[counterI+1] % 6;
	    SPoints[B][number]->get_point_ud(dir)->set_potential_point(dataR_in[counterI/2]);
	    counterI+=2;
	  }
      }
  }
}
