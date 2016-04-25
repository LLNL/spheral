#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void add_buffer_values(Fractal_Memory& mem,int level,
			 vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    static int _COUNTER=0;
    int FractalRank=mem.p_mess->FractalRank;
    bool Ranky=FractalRank==33;
    cerr << " SOLVED D " << _COUNTER << " " << FractalRank << "\n";
    vector <int>BOX=mem.BoxesLev[FractalRank][level];
    // vector <int>BBOX=mem.BBoxesLev[FractalRank][level];
    int HypreNodes=mem.p_mess->HypreNodes;
    vector <int> counts_in(HypreNodes);
    vector <int> counts_out(HypreNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(HypreNodes);
    vector < vector <double> > dataR_out(HypreNodes);
    vector <int>pos(3);
    for(int FR : mem.Touchy)
      {
	int HR=mem.p_mess->IHranks[FR];
	assert(HR >= 0);
	vector <int>FRPBOX=mem.PBoxesLev[FR][level];
	vector <int>FRBOX=mem.BoxesLev[FR][level];
	int B=0;
	for(vector <int>& SB : SBoxes)
	  {
	    if(overlap_boxes(FRPBOX,SB))
	      {
		for(Point* &p : SPoints[B])
		  {
		    p->get_pos_point(pos);
		    if(vector_in_box(pos,FRPBOX) && !vector_in_box(pos,FRBOX))
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
	cerr << " OUTADD " << _COUNTER << " " << FractalRank << " " << FR << " " << HR << " " << counts_out[HR] << "\n";
      }
    cerr << " SOLVED E " << _COUNTER << " " << FractalRank << "\n";
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=3;
    int doubles=1;
    mem.p_mess->Send_Data_Some_How(9,mem.p_mess->HypreWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    cerr << " SOLVED F " << _COUNTER << " " << FractalRank << " " << BOX[0]  << " " << BOX[1]  << " " << BOX[2]  << " " << BOX[3]  << " " << BOX[4]  << " " << BOX[5] << "\n";
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    array<int,3>ar3;
    std::map<array<int,3>,Point*,point_comp> edgeP;
    // std::pair<std::map<array<int,3>,Point*>::iterator,bool> succ;
    int inP=0;
    for(Group* &pg : mem.all_groups[level])
      {
	if(pg->get_buffer_group())
	  {
	    for(Point* &p : pg->list_points)
	      {
		p->get_pos_point(pos);
		if(!vector_in_box(pos,BOX))
		  {
		    std::move(pos.begin(),pos.end(),ar3.begin());
		    // succ=edgeP.insert(pap(ar3,p));
		    edgeP[ar3]=p;
		    inP++;
		  }
	      }
	  }
      }
    cerr << " SOLVED G " << _COUNTER << " " << FractalRank << " " << inP << " " << edgeP.size() << "\n";
    int c1=0;
    int c3=0;
    int found=0;
    int notfound=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	for(int c=0;c<counts_in[HR];c++)
	  {
	    array <int,3>posin={dataI_in[c3],dataI_in[c3+1],dataI_in[c3+2]};
	    if(Ranky)
	      cerr << " EdgePosB " << FractalRank << " " << level << " " << dataI_in[c3]  << " " << dataI_in[c3+1]  << " " << dataI_in[c3+2] << "\n";
	    // std::map<char,int>::iterator it;
	    std::map<array<int,3>,Point*>::iterator eb=edgeP.find(posin);
	    if(eb != edgeP.end())
	      {
		// Point* PF=eb->second;
		// PF->set_potential_point(dataR_in[c1]);
		eb->second->set_potential_point(dataR_in[c1]);
		found++;
	      }
	    else
	      {
		notfound++;
	      }
	    c1++;
	    c3+=3;
	  }
      }
    cerr << " SOLVED H " << _COUNTER << " " << FractalRank << "\n";
    cerr << " ADDUP " << _COUNTER << " " << FractalRank << " " << found << " " << notfound << "\n";
    _COUNTER++;
  }
}
