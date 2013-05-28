#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool hypre_ij_numbering(Fractal_Memory& mem,Fractal& frac,vector <Point*>& hypre_points,const int& level)
  {
    frac.timing(-1,32);
    const int FractalRank=mem.p_mess->FractalRank;
    const int FractalNodes=mem.p_mess->FractalNodes;
    ofstream& FH=mem.p_file->FileHypre;
    const int zoom=Misc::pow(2,frac.get_level_max());
    const int spacing=Misc::pow(2,frac.get_level_max()-level);
    const int mult=Misc::pow(2,level);
    int wrap=zoom*frac.get_grid_length();
    const bool period=frac.get_periodic();
    vector <int> posa(3);
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	if(!group.get_buffer_group() && group.list_points.size() < mem.min_hypre_group_size)
	  continue;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  hypre_points.push_back(*point_itr);
      }
    int count=hypre_points.size();
    mem.p_mess->IAmAHypreNode=count > 0;
    /*
    if(count==0)
      {
	count=1;
	hypre_points.push_back(0);
	FH << " fake a" << endl;
      }
    */
    frac.timing(-1,36);
    mem.p_mess->Full_Stop();
    frac.timing(1,36);
    mem.p_mess->How_Many_On_Nodes(count,mem.ij_counts);
    mem.ij_offsets.resize(FractalNodes);
    mem.ij_offsets[0]=0;

    vector <int> ranks(FractalNodes);
    int nhyp=0;
    for(int ni=0;ni<FractalNodes;ni++)
      {
	if(mem.ij_counts[ni] == 0)
	  continue;
	ranks[nhyp]=ni;
	if(ni > 0)
	  mem.ij_offsets[nhyp]=mem.ij_offsets[ni-1]+mem.ij_counts[ni-1];
	nhyp++;
      }
    int HypreNodes=nhyp;
    if(HypreNodes == 0)
      {
	frac.timing(1,32);
	return false;
      }
    mem.p_mess->HypreNodes=HypreNodes;
    for(int ni=0;ni<HypreNodes;ni++)
      FH << " offsets " << ni << " " << level << " " << mem.ij_offsets[ni] << " " << mem.ij_counts[ni] << endl;
    int totals=mem.ij_offsets[FractalNodes-1]+mem.ij_counts[FractalNodes-1];
    mem.p_mess->HypreGroupCreate(ranks);
    int HypreRank=mem.p_mess->what_is_my_Hypre_rank();
    mem.p_mess->HypreRank=HypreRank;
    assert(HypreNodes == mem.p_mess->how_many_Hypre_nodes());
    FH << " total " << totals << " " << level << endl;
    vector <int>PBox(6);
    vector <int>PBoxLeft(3);
    PBox=mem.PBoxesLev[FractalRank][level];
    PBoxLeft[0]=PBox[0];
    PBoxLeft[1]=PBox[2];
    PBoxLeft[2]=PBox[4];
    int lengthx=PBox[1]-PBox[0]+1;
    int lengthy=PBox[3]-PBox[2]+1;
    int lengthz=PBox[5]-PBox[4]+1;
    vector < vector < vector <Point*> > > searches;
    searches.resize(3);
    searches[0].resize(lengthx);
    searches[1].resize(lengthy);
    searches[2].resize(lengthz);
    count=mem.ij_offsets[FractalRank];
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point*p=*point_itr;
	if(p)p->set_ij_number(count);
	count++;
      }
    vector <int>pp;
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point*p=*point_itr;
	if(p == 0) continue;
	p->set_ij_neighbors();
      }
    vector <Point*> send_list;
    int x,y,z;
    vector <int>pos(3);
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p == 0) continue;
	if(p->to_be_sent())
	  {
	    send_list.push_back(p);
	  }
	else if(p->to_receive())
	  {
	    p->get_pos_point(pos);
	    x=(pos[0]-PBoxLeft[0])/spacing;
	    assert(x >= 0 && x < lengthx);
	    searches[0][x].push_back(p);
	    y=(pos[1]-PBoxLeft[1])/spacing;
	    assert(y >= 0 && y < lengthy);
	    searches[1][y].push_back(p);
	    z=(pos[2]-PBoxLeft[2])/spacing;
	    assert(z >= 0 && z < lengthz);
	    searches[2][z].push_back(p);
	  }
      }
    vector <int> pos_lefts(3);
    vector <int> pos_rights(3);
    if(period)
      left_right(send_list,pos_lefts,pos_rights,wrap);
    else
      left_right(send_list,pos_lefts,pos_rights);
    vector <int> counts_out(FractalNodes);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    int ssize=send_list.size();
    for(int FR=0;FR<FractalNodes;FR++)
      {
	counts_out[FR]=0;
	if(FR == FractalRank)
	  continue;
	PBox=mem.PBoxesLev[FR][level];
	if(!overlap(pos_lefts,pos_rights,PBox))
	  continue;
	for(int ni=0;ni<ssize;ni++)
	  {
	    send_list[ni]->get_pos_point(pos);
	    if(period)
	      {
		for(int ni=0;ni<3;ni++)
		  pos[ni]=(pos[ni]+wrap) % wrap;
	      }
	    if(!overlap(pos,pos,PBox))
	      continue;
	    dataI_out[FR].push_back(pos[0]);
	    dataI_out[FR].push_back(pos[1]);
	    dataI_out[FR].push_back(pos[2]);
	    dataI_out[FR].push_back(send_list[ni]->get_ij_number());
	    counts_out[FR]++;
	  }
      }
    vector <int> counts_in(FractalNodes);
    vector <double> dataR_in;
    vector <int> dataI_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=4;
    int doubles=0;
    frac.timing(-1,37);
    mem.p_mess->Full_Stop();
    frac.timing(1,37);
    mem.p_file->note(true," hypre numbering a ");
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    mem.p_file->note(true," hypre numbering b ");
    mem.p_mess->Send_Data_Somewhere(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," hypre numbering c ");
    dataR_out.clear();
    dataI_out.clear();
    Point* p=0;
    vector <int>xyz(3);
    int ni4=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    xyz[0]=(dataI_in[ni4  ]-PBoxLeft[0])/spacing;
	    xyz[1]=(dataI_in[ni4+1]-PBoxLeft[1])/spacing;
	    xyz[2]=(dataI_in[ni4+2]-PBoxLeft[2])/spacing;
	    int ni=shortest_vector(searches[0][xyz[0]],searches[1][xyz[1]],searches[2][xyz[2]]);
	    int number=which_element(searches[ni][xyz[ni]],dataI_in[ni4],dataI_in[ni4+1],dataI_in[ni4+2],FH);
	    if(number < 0)
	      FH << "not found " << FR << " " << c << endl;
	    else
	      searches[ni][xyz[ni]][number]->copy_ij_index(dataI_in[ni4+3]);
	    ni4+=4;
	  }
      }
    frac.timing(1,32);
    return hypre_points.size() > 0;
  }
}
