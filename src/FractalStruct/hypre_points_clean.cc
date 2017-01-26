#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_clean(Fractal_Memory& mem,int level,vector< vector<Point*> >& hypre_points)
  {
    // ofstream& FHT=mem.p_file->DUMPS;
    vector <int>Touchies;
    for(int FR : mem.TouchWhichBoxes)
      if(mem.p_mess->counts_on_nodes[level][2*FR+1])
	{
	  // FHT << " Touchies " << level << " " << FR << "\n";
	  Touchies.push_back(FR);
	}

    int FractalRank=mem.p_mess->FractalRank;      
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <int> counts_in(FractalNodes);
    vector <int> counts_out(FractalNodes,0);
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector < vector <int> > dataI_out(FractalNodes);
    vector < vector <double> > dataR_out(FractalNodes);

    vector<Point*> edgies;
    vector <int>BOX=mem.BoxesLev[FractalRank][level];
    vector <int>pos(3);
    vector <int>pos1(3);
    array <int,3>ar3;
    for(auto pg : hypre_points)
      for(auto p : pg)
	{
	  if(on_edge(p,BOX))
	    {
	      edgies.push_back(p);
	      for(int ni=0;ni<6;ni++)
		{
		  Point* p1=p->get_point_ud_0(ni,-13);
		  if(vector_in_box(p1,BOX))
		    continue;
		  p1->get_pos_point(pos1);
		  for(int FR : Touchies)
		    {
		      if(on_edge(pos1,mem.BoxesLev[FR][level]))
			{
			  dataI_out[FR].push_back(pos1[0]);
			  dataI_out[FR].push_back(pos1[1]);
			  dataI_out[FR].push_back(pos1[2]);
			  dataI_out[FR].push_back(ni);
			  counts_out[FR]++;
			  break;
			}
		    }
		}
	    }
	}
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=4;
    int doubles=0;
    mem.p_mess->Send_Data_Some_How(12,mem.p_mess->FractalWorld,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    counts_out.clear();
    dataI_out.clear();
    dataR_out.clear();
    dataR_in.clear();
    std::multimap<array<int,3>,int,point_comp2> edgein;
    int c4=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    array <int,3>posin={dataI_in[c4],dataI_in[c4+1],dataI_in[c4+2]};
	    assert(on_edge(posin,BOX));
	    // FHT << " POSIN " << posin[0] << " "  << posin[1] << " "  << posin[2] << "\n";
	    edgein.insert(pair<array<int,3>,int>(posin,dataI_in[c4+3]));
	    c4+=4;
	  }
      }
    dataI_in.clear();
    for(auto &p : edgies)
      {
	p->get_pos_point(ar3);
	int count=edgein.count(ar3);
	if(count == 0)
	  continue;
	auto found=edgein.equal_range(ar3);
	// FHT << " FOUND " << count << " " << ar3[0] << " " << ar3[1] << " " << ar3[2] << " ";
	for(int ni=0;ni<6;ni++)
	  {
	    Point* p1=p->get_point_ud_0(ni,18);
	    if(vector_in_box(p1,BOX))
	      continue;
	    int ni1=(ni % 2) == 1 ? ni-1:ni+1;
	    bool success=false;
	    for(auto f=found.first;f != found.second;f++)
	      {
		success=f->second == ni1;
		if(success)
		  break;
	      }
	    if(!success)
	      p->set_trouble(true);
	    // FHT << ni << " " << ni1 << " " << success << p->get_trouble() << "\n";
	  }
      }
    
  }
}
