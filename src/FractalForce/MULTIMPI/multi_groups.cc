#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void multi_groups(Fractal_Memory& mem,int& levelmin,int& levelmax)
  {
    static bool donezero=false;
    vector <int>pos(3);
    int TWB=mem.TouchWhichBoxes.size();
    if(levelmin == 0)
      {
	in zoom=Misc::pow(2,mem.level_max);
	Group* p_group=mem.all_groups[0][0];
	Group& group=*p_group;
	for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  {
	    Point* p_point=*point_itr;
	    p_point->clear_ghosts();
	    if(!p_point->get_buffer_point() && !p_point->get_edge_point())
	      continue;
	    p_point->get_pos_point(pos);
	    pos[0]/=zoom;
	    pos[1]/=zoom;
	    pos[2]/=zoom;
	    for (int TW=0;TW<TWB;TW++)
	      {
		int FR=TouchWhichBoxes[TW];
		if(!vector_in_box(pos,mem.BBoxes[FR]))
		  continue;
		int nx=pos[0]-mem.PBoxes[FR][0];
		int ny=pos[1]-mem.PBoxes[FR][1];
		int nz=pos[2]-mem.PBoxes[FR][2];
		in n=nx+(ny+nz*PBoxesLength[2])*PBoxesLength[1];
		p_point->add_ghost(FR,0,n);
	      }
	  }
	donezero=true;
	if(levelmax == 0)
	  return;
      }
    assert(donezero);
    vector <int>ghosts_FR;
    vector <int>ghosts_group_number;
    vector <int>ghosts_point_number;
    int rp=-1;
    vector <int> counts_out;
    vector <vector <int> > dataI_out;
    vector <vector <double> > dataR_out;
    for (int level=1;level <= levelmax;level++)
      {
	counts_out.assign(FractalNodes,0);
	dataI_out.clear();
	dataI_out.resize(FractalNodes);
	dataR_out.clear();
	dataR_out.resize(FractalNodes);
	for(vector <Group*>::const_iterator group_itr=_mem.all_buffer_groups[level].begin();
	    group_itr!=mem.all_buffer_groups[level].end();group_itr++)
	  {
	    Group* p_group=*group_itr;
	    Group& group=*p_group;
	    int group_number=p_group->get_number_group();
	    int point_number=-1;
	    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	      {
		Point* p_point=*point_itr;
		p_point->clear_ghosts();
		point_number++;
		Point* p_Point=p_point->move_rp(0)->get_point_pointer();
		if(!p_Point->any_ghosts())
		  continue;
		p_point->get_real_pointer(rp);
		p_Point->get_ghosts(ghosts_FR,ghosts_group_number,ghosts_point_number);
		int FRNodes=ghosts_FR.size();
		for(int node=0;node<FRNodes;node++)
		  {
		    int FR=ghosts_FR[node];
		    dataI_out[FR].push_back(ghosts_group_number[node]);
		    dataI_out[FR].push_back(ghosts_point_number[node]);
		    dataI_out[FR].push_back(group_number);
		    dataI_out[FR].push_back(point_number);
		    dataI_out[FR].push_back(rp);
		    counts_out[FR]++;
		  }
	      }
	  }
	vector <int> counts_in(FractalNodes);
	vector <double> dataR_in;
	vector <int> dataI_in;
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=5;
	int doubles=0;
	mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
	mem.p_mess->Send_Data_Somewhere_No_Block(counts_out,counts_in,integers,doubles,
						 dataI_out,dataI_in,how_manyI,
						 dataR_out,dataR_in,how_manyR);
	dataR_out.clear();
	dataI_out.clear();
	int ni5=0;
	for(int TW=0;TW<TWB;TW++)
	  {
	    int FR=mem.TouchWhichBoxes[TW];
	    for(int c=0;c<counts_in[FR];c++)
	      {
		Group* p_group=mem.all_groups[level-1][dataI_in[ni5]];
		Point* p_point=p_group->list_points[dataI_in[ni5+1]];
		p_point=p_point->get_daughter_point();
		p_point=p_point->move_rp(dataI_in[ni5+4]);
		p_point->set_ghosts(FR,dataI_in[ni5+2],dataI_in[ni5+3]);
		//mem.all_groups[level-1][dataI_in[ni5]]->list_points[dataI_in[ni5+1]]->get_daughter_point()->move_rp(dataI_in[ni5+4])->set_ghosts(FR,dataI_in[ni5+2],dataI_in[ni5+3]);
		ni5+=5;
	      }
	  }
	counts_in.assign(FractalNodes,0);
	dataI_in.clear();
	dataR_in.clear();
      }
    int lmin=max(1,levelmin);
    for(int level=lmin,level<=levelmax)
      {
	for(vector <Group*>::const_iterator group_itr=_mem.all_buffer_groups[level].begin();
	    group_itr!=mem.all_buffer_groups[level].end();group_itr++)
	  {
	    Group* p_group=*group_itr;
	    Group& group=*p_group;
	    int group_number=p_group->get_number_group();
	    int point_number=-1;
	    for( vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	      {
		Point* p_point=*point_itr;
	      }	
      }
  }
}
