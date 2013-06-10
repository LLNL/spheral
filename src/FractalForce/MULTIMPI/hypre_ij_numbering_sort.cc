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
    bool period=frac.get_periodic();
    vector <int> posa(3);
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  {
	    Point* p=*point_itr;
	    if(!p->get_really_passive())
	      hypre_points.push_back(p);
	    //	    else
	    //	      FH << "really " << p->get_pos_point(0) << " " << p->get_pos_point(1) << " " << p->get_pos_point(2) << endl;
	  }    
      }
    int count=hypre_points.size();
    if(count==0)
      {
	count=1;
	hypre_points.push_back(0);
	FH << " fake a" << endl;
      }
    frac.timing(-1,36);
    mem.p_mess->Full_Stop();
    frac.timing(1,36);
    mem.p_mess->How_Many_On_Nodes(count,mem.ij_counts);
    mem.ij_offsets.resize(FractalNodes);
    mem.ij_offsets[0]=0;
    for(int ni=1;ni<FractalNodes;ni++)
      mem.ij_offsets[ni]=mem.ij_offsets[ni-1]+mem.ij_counts[ni-1];
    for(int ni=0;ni<FractalNodes;ni++)
      FH << " offsets " << ni << " " << level << " " << mem.ij_offsets[ni] << " " << mem.ij_counts[ni] << endl;
    int totals=mem.ij_offsets[FractalNodes-1]+mem.ij_counts[FractalNodes-1];
    if(totals <= FractalNodes)
      {
	frac.timing(1,32);
	return false;
      }
    FH << " total " << totals << " " << level << endl;
    vector <int>PBox(6);
    count=mem.ij_offsets[FractalRank];
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point*p=*point_itr;
	if(p)
	  p->set_ij_number(count);
	count++;
      }
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p) 
	  p->set_ij_neighbors();
      }
    vector <Point*> send_list;
    vector <Point*> receive_list;
    vector <int>pos(3);
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p == 0) 
	  continue;
	if(p->to_be_sent())
	  send_list.push_back(p);
	else if(p->to_receive())
	  receive_list.push_back(p);
      }
    sort3_list(receive_list,0);
    /*
    for(vector<Point*>::const_iterator point_itr=receive_list.begin();point_itr !=receive_list.end();++point_itr)
      {
	Point* p=*point_itr;
	FH << "receive " << p->get_pos_point(0) << " " << p->get_pos_point(1) << " " << p->get_pos_point(2) << endl;
      }
    */
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
	    //	    FH << " send " << FR << " " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
	    dataI_out[FR].push_back(pos[0]);
	    dataI_out[FR].push_back(pos[1]);
	    dataI_out[FR].push_back(pos[2]);
	    dataI_out[FR].push_back(send_list[ni]->get_ij_number());
	    counts_out[FR]++;
	  }
      }
    send_list.clear();
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
    Point* psend= new Point;
    int what=0;
    Misc::dim2=(what+2) % 3;
    Misc::dim1=(what+1) % 3;
    Misc::dim0=(what+0) % 3;
    mem.p_mess->Full_Stop();
    int ni4=0;
    int found=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    //	    FH << " try this " << FR << " " << c << " " << ni4 << endl;
	    //	    FH << dataI_in[ni4] << " " << dataI_in[ni4+1] << " " << dataI_in[ni4+2] << " " << psend << endl;
	    psend->set_pos_point(dataI_in[ni4],dataI_in[ni4+1],dataI_in[ni4+2]);
	    int number=findPointByPos(receive_list,psend,FH);
	    if(number < 0)
	      FH << "not found OK " << FR << " " << c << " " << dataI_in[ni4] << " " << dataI_in[ni4+1] << " " << dataI_in[ni4+2] << endl;
	    else
	    if(number >= 0)
	      {
		receive_list[number]->copy_ij_index(dataI_in[ni4+3]);
		found++;
		//		FH << "found OK " << FR << " " << c << " " << dataI_in[ni4] << " " << dataI_in[ni4+1] << " " << dataI_in[ni4+2] << endl;
	      }
	    ni4+=4;
	  }
      }
    delete psend;
    FH << "found connections " << found << endl;
    /*
    for(vector<Point*>::const_iterator point_itr=receive_list.begin();point_itr !=receive_list.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p->get_found_it())
	  continue;
	FH << " not found really bad" << endl;
	p->dumpp(FH);
      }
    */
    //    mem.p_mess->Full_Stop();
    //    assert(0);
    frac.timing(1,32);
    return true;
  }
}
