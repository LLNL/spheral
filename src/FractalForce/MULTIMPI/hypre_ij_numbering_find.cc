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
    vector <int>pos(3);
    vector <int>HBox=mem.HRBoxesLev[FractalRank][level];
    vector <int>HRBox=HBox;
    vector <int>HSBox=mem.HSBoxesLev[FractalRank][level];
    FH << " HBox " << HBox[0] << " "  << HBox[1] << " "  << HBox[2] << " "  << HBox[3] << " "  << HBox[4] << " "  << HBox[5] << endl;
    FH << " HSBox " << HSBox[0] << " "  << HSBox[1] << " "  << HSBox[2] << " "  << HSBox[3] << " "  << HSBox[4] << " "  << HSBox[5] << endl;
    vector <Point*>hypre_pointsB;
    unsigned int minsize=mem.min_hypre_group_size;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	if(group.list_points.size() <= minsize && !group.get_buffer_group())
	  continue;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  hypre_points.push_back(*point_itr);
	if(!group.get_buffer_group())
	  continue;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  hypre_pointsB.push_back(*point_itr);
      }
    int count=hypre_points.size();
    if(count==0)
      {
	count=1;
	hypre_points.push_back(0);
	FH << " fake a" << endl;
      }
    mem.p_mess->How_Many_On_Nodes(count,mem.ij_counts);
    mem.ij_counts.resize(FractalNodes+1);
    mem.ij_offsets.assign(FractalNodes+1,0);
    FH << " offsets numbering " << 0 << " " << level << " " << mem.ij_offsets[0] << " " << mem.ij_counts[0] << endl;
    for(int FR=1;FR<=FractalNodes;FR++)
      {
	mem.ij_offsets[FR]=mem.ij_offsets[FR-1]+mem.ij_counts[FR-1];
	FH << " offsets numbering " << FR << " " << level << " " << mem.ij_offsets[FR] << " " << mem.ij_counts[FR] << endl;
      }
    mem.ij_countsB=mem.ij_counts;
    mem.ij_offsetsB=mem.ij_offsets;
    int totals=mem.ij_offsets[FractalNodes];
    assert(totals >= FractalNodes);
    if(totals == FractalNodes)
      {
	frac.timing(1,32);
	return false;
      }
    //    FH << " total points " << totals << " " << involved << " " << nodes_eff << " " << loading << " " << count_max << " " << level << endl;
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
	  p->set_ij_neighbors(HRBox);
      }
    vector <Point*> send_list;
    vector <Point*> receive_list;
    for(vector<Point*>::const_iterator point_itr=hypre_pointsB.begin();point_itr !=hypre_pointsB.end();++point_itr)
      {
	Point* p=*point_itr;
	p->get_pos_point(pos);
	if(on_edge(pos,HSBox))
	  send_list.push_back(p);
	if(on_edge(pos,HRBox))
	  receive_list.push_back(p);
      }
    sort3_list(receive_list,0);
    vector <int> pos_lefts(3);
    vector <int> pos_rights(3);
    left_right(send_list,pos_lefts,pos_rights);
    vector <int> counts_out(FractalNodes,0);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    int ssize=send_list.size();
    int TWB=mem.TouchWhichBoxes.size();
    for(int TW=0;TW<TWB;TW++)
      {
	int FR=mem.TouchWhichBoxes[TW];
	counts_out[FR]=0;
	if(FR == FractalRank)
	  continue;
	HRBox=mem.HRBoxesLev[FR][level];
	if(!overlap(pos_lefts,pos_rights,HRBox))
	  continue;
	for(int ni=0;ni<ssize;ni++)
	  {
	    send_list[ni]->get_pos_point(pos);
	    if(!vector_in_box(pos,HRBox))
	      continue;
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
    mem.p_file->note(true," hypre numbering a ");
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    mem.p_file->note(true," hypre numbering b ");
    mem.p_mess->Send_Data_Somewhere_No_Block(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," hypre numbering c ");
    dataR_out.clear();
    dataI_out.clear();
    Point* psend= new Point;
    int what=0;
    Misc::dim2=(what+2) % 3;
    Misc::dim1=(what+1) % 3;
    Misc::dim0=(what+0) % 3;
    int ni4=0;
    int found=0;
    for(int TW=0;TW<TWB;TW++)
      {
	int FR=mem.TouchWhichBoxes[TW];
	for(int c=0;c<counts_in[FR];c++)
	  {
	    psend->set_pos_point(dataI_in[ni4],dataI_in[ni4+1],dataI_in[ni4+2]);
	    int number=findPointByPos(receive_list,psend,FH);
	    if(number < 0)
	      FH << "not found OK " << FR << " " << c << " " << dataI_in[ni4] << " " << dataI_in[ni4+1] << " " << dataI_in[ni4+2] << endl;
	    else
	      {
		receive_list[number]->copy_ij_index(dataI_in[ni4+3]);
		found++;
	      }
	    ni4+=4;
	  }
      }
    delete psend;
    FH << "found connections " << FractalRank << " " << level << " " << found << endl;
    frac.timing(1,32);
    return true;
  }
}
