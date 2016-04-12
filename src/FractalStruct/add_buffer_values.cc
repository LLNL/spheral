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
    bool Ranky=FractalRank==32;
    cerr << " SOLVED D " << _COUNTER << " " << FractalRank << endl;
    vector <int>BOX=mem.BoxesLev[FractalRank][level];
    vector <int>BBOX=mem.BBoxesLev[FractalRank][level];
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
	vector <int>PBOX=mem.BBoxesLev[FR][level];
	int B=0;
	for(vector <int>& SB : SBoxes)
	  {
	    if(overlap_boxes(PBOX,SB))
	      {
		for(Point* &p : SPoints[B])
		  {
		    p->get_pos_point(pos);
		    if(vector_in_box(pos,PBOX))
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
	cerr << " OUTADD " << _COUNTER << " " << FractalRank << " " << FR << " " << HR << " " << counts_out[HR] << endl;
      }
    cerr << " SOLVED E " << _COUNTER << " " << FractalRank << endl;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=3;
    int doubles=1;
    mem.p_mess->Send_Data_Some_How(9,mem.p_mess->HypreWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    cerr << " SOLVED F " << _COUNTER << " " << FractalRank << " " << BOX[0]  << " " << BOX[1]  << " " << BOX[2]  << " " << BOX[3]  << " " << BOX[4]  << " " << BOX[5] << endl;
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    std::set<Point*>edgeP;
    std::pair<std::set<Point*>::iterator,bool> succ;
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
		    succ=edgeP.insert(p);
		    assert(succ.second);
		    // if(Ranky)
		    //   cerr << " EDGEA " << _COUNTER << " " << FractalRank << " " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
		    inP++;
		  }
	      }
	  }
      }
    cerr << " SOLVED G " << _COUNTER << " " << FractalRank << " " << inP << " " << edgeP.size() << endl;
    auto iend=edgeP.end();
    int c1=0;
    int c3=0;
    int found=0;
    int notfound=0;
    Point* p=new Point;
    vector <int>posin(3);
    for(int HR=0;HR<HypreNodes;HR++)
      {
	for(int c=0;c<counts_in[HR];c++)
	  {
	    p->set_pos_point(dataI_in[c3],dataI_in[c3+1],dataI_in[c3+2]);
	    p->get_pos_point(posin);
	    // assert(vector_in_box(posin,BBOX));
	    // assert(!vector_in_box(posin,BOX));
	    auto it=edgeP.find(p);
	    if(it != iend)
	      {
		if(Ranky)
		  cerr << " ITISFOUND " ;
		(*it)->set_potential_point(dataR_in[c1]);
		found++;
	      }
	    else
	      {
		if(Ranky)
		  cerr << " NOTFOUND " ;
		notfound++;
	      }
	    if(Ranky)
	      cerr << _COUNTER << " " << FractalRank << " " << dataI_in[c3] << " " << dataI_in[c3+1] << " " << dataI_in[c3+2] << endl;
	    c1++;
	    c3+=3;
	  }
      }
    delete p;
    cerr << " SOLVED H " << _COUNTER << " " << FractalRank << endl;
    cerr << " ADDUP " << _COUNTER << " " << FractalRank << " " << found << " " << notfound << "\n";
    _COUNTER++;
  }
}
