#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool hypre_ij_numbering(Fractal_Memory& mem,HypHest& HYP,vector <Point*>& hypre_points,const int& level)
  {
    double time0=mem.p_mess->Clock();
    const int FractalRank=mem.p_mess->FractalRank;
    const int FractalNodes=mem.p_mess->FractalNodes;
    FILE* PFH=mem.p_file->PFHypre;
    vector <int>pos(3);
    vector <int>HBox=mem.HRBoxesLev[FractalRank][level];
    vector <int>HRBox=HBox;
    vector <int>HSBox=mem.HSBoxesLev[FractalRank][level];
    vector <Point*> send_list;
    vector <Point*> receive_list;
    fprintf(PFH," HBox %d %d %d %d %d %d \n",HBox[0],HBox[1],HBox[2],HBox[3],HBox[4],HBox[5]);
    fprintf(PFH," HSBox %d %d %d %d %d %d \n",HSBox[0],HSBox[1],HSBox[2],HSBox[3],HSBox[4],HSBox[5]);
    //    vector <Point*>hypre_pointsB;
    unsigned int minsize=mem.min_hypre_group_size;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	bool bg=group.get_buffer_group();
	if(group.list_points.size() <= minsize && !bg)
	  continue;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  {
	    Point* p=*point_itr;
	    p->get_pos_point(pos);
	    if(p->get_inside() || on_edge(pos,HSBox) || on_edge(pos,HRBox))
	      hypre_points.push_back(p);
	  }
	if(!bg)
	  continue;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  {
	    Point* p=*point_itr;
	    p->get_pos_point(pos);
	    if(on_edge(pos,HSBox))
	      send_list.push_back(p);
	    if(on_edge(pos,HRBox))
	      receive_list.push_back(p);
	  }
      }
    double time1=mem.p_mess->Clock();
    int count=hypre_points.size();
    mem.p_mess->IAmAHypreNode=count > 0;
//     mem.p_mess->Full_Stop_Do_Not_Argue();
//     double tguess=-mem.p_mess->Clock();
    mem.p_mess->How_Many_On_Nodes(count,HYP.ij_counts);
//     tguess+=mem.p_mess->Clock();
    mem.Touchy=mem.TouchWhichBoxes;
    //
    vector <int> pos_lefts(3);
    vector <int> pos_rights(3);
    left_right(send_list,pos_lefts,pos_rights);
    mem.Touchy.clear();
    int TWB=mem.TouchWhichBoxes.size();
    for(int TW=0;TW<TWB;TW++)
      {
	if(!mem.p_mess->IAmAHypreNode)
	  continue;
	int FR=mem.TouchWhichBoxes[TW];
	if(FR != FractalRank && HYP.ij_counts[FR] > 0 && overlap(pos_lefts,pos_rights,mem.HRBoxesLev[FR][level]))
	  mem.Touchy.push_back(FR);
      }    
    int totals=std::accumulate(HYP.ij_counts.begin(),HYP.ij_counts.end(),0);
    double time2=mem.p_mess->Clock();
    node_groups(mem,HYP);
    HYP.ij_counts.resize(FractalNodes);
    HYP.ij_offsets.assign(FractalNodes+1,0);
    mem.p_mess->Hranks.clear();
    mem.p_mess->IHranks.assign(FractalNodes,-1);
    int HypreNodes=0;
    for(int FR=0;FR<=FractalNodes;FR++)
      {
	if(FR > 0)
	  HYP.ij_offsets[FR]=HYP.ij_offsets[FR-1]+HYP.ij_counts[FR-1];
	if(FR < FractalNodes && HYP.ij_counts[FR] > 0)
	  {
	    mem.p_mess->Hranks.push_back(FR);
	    mem.p_mess->IHranks[FR]=HypreNodes;
	    HypreNodes++;
	  }
      }
    if(totals == 0)
      return false;
    int HR=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	if(HYP.ij_counts[FR] == 0)
	  continue;
	HYP.ij_counts[HR]=HYP.ij_counts[FR];
	HYP.ij_offsets[HR]=HYP.ij_offsets[FR];
	HR++;
      }
    double time3=mem.p_mess->Clock();
    mem.p_mess->HypreNodes=HypreNodes;
    mem.p_mess->HypreGroupCreate(mem.p_mess->Hranks);
    int HypreRank=mem.p_mess->MyHypreRank();
    if(mem.p_mess->IAmAHypreNode && FractalRank == mem.Touchy.front())
      {
	for(int TW=0;TW<mem.Touchy.size();TW++)
 	  {
	    // 	    cerr << " HRANKS " << mem.steps << " " << level << " " << FractalRank << " " << TW << " " << mem.Touchy[TW] << "\n";
 	    fprintf(PFH," HRANKS %d %d %d %d %d \n",mem.steps,level,FractalRank,TW,mem.Touchy[TW]);
 	  }
      }
    HYP.ij_counts.resize(HypreNodes);
    HYP.ij_offsets.resize(HypreNodes+1);
    if(mem.p_mess->IAmAHypreNode)
      HYP.ij_offsets[HypreNodes]=HYP.ij_offsets[HypreNodes-1]+HYP.ij_counts[HypreNodes-1];
    mem.p_mess->HypreRank=mem.p_mess->MyHypreRank();
    assert(HypreRank == mem.p_mess->HypreRank);
    HYP.ij_countsB=HYP.ij_counts;
    HYP.ij_offsetsB=HYP.ij_offsets;
    if(!mem.p_mess->IAmAHypreNode)
      return true;
    count=HYP.ij_offsets[HypreRank];
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point*p=*point_itr;
	p->set_ij_number(count);
	count++;
      }
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	Point* p=*point_itr;
	p->set_ij_neighbors(HRBox);
      }
    sort3_list(receive_list,0);
    double time4=mem.p_mess->Clock();
    vector <int> counts_out(HypreNodes,0);
    vector <vector <int> > dataI_out(HypreNodes);
    vector <vector <double> > dataR_out(HypreNodes);
    int ssize=send_list.size();
    TWB=mem.TouchWhichBoxes.size();
    for(int TW=0;TW<TWB;TW++)
      {
	int FR=mem.TouchWhichBoxes[TW];
	if(FR == FractalRank)
	  continue;
	int HR=mem.p_mess->IHranks[FR];
	if(HR < 0)
	  continue;
	counts_out[HR]=0;
	HRBox=mem.HRBoxesLev[FR][level];
	if(!overlap(pos_lefts,pos_rights,HRBox))
	  continue;
	for(int ni=0;ni<ssize;ni++)
	  {
	    send_list[ni]->get_pos_point(pos);
	    if(!vector_in_box(pos,HRBox))
	      continue;
	    dataI_out[HR].push_back(pos[0]);
	    dataI_out[HR].push_back(pos[1]);
	    dataI_out[HR].push_back(pos[2]);
	    dataI_out[HR].push_back(send_list[ni]->get_ij_number());
	    counts_out[HR]++;
	  }
      }
    send_list.clear();
    double time5=mem.p_mess->Clock();
    vector <int> counts_in(HypreNodes);
    vector <double> dataR_in;
    vector <int> dataI_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=4;
    int doubles=0;
    double time6=mem.p_mess->Clock();
    mem.p_file->note(true," hypre numbering b ");
    mem.p_mess->Send_Data_Some_How(2,mem.p_mess->HypreWorld,
				   counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    double time7=mem.p_mess->Clock();
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
    TWB=mem.TouchWhichBoxes.size();
    for(int TW=0;TW<TWB;TW++)
      {
	int FR=mem.TouchWhichBoxes[TW];
	int HR=mem.p_mess->IHranks[FR];
	if(HR < 0)
	  continue;
	for(int c=0;c<counts_in[HR];c++)
	  {
	    psend->set_pos_point(dataI_in[ni4],dataI_in[ni4+1],dataI_in[ni4+2]);
	    int number=findPointByPos(receive_list,psend,PFH);
	    if(number < 0)
	      fprintf(PFH,"not found OK %d %d %d %d %d \n",FR,c,dataI_in[ni4],dataI_in[ni4+1],dataI_in[ni4+2]);
	    else
	      {
		receive_list[number]->copy_ij_index(dataI_in[ni4+3]);
		found++;
	      }
	    ni4+=4;
	  }
      }
    double time8=mem.p_mess->Clock();
    delete psend;
    fprintf(PFH,"found connections %d %d %d \n",FractalRank,level,found);
    fprintf(mem.p_file->PFTime," Hypre Find %3d %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E \n",
	    level,time1-time0,time2-time1,time3-time2,time4-time3,time5-time4,time6-time5,time7-time6,time8-time7,time8-time0);
    return true;
  }
}
