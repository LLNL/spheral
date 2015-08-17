#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool hypre_ij_numbering_selfie(Fractal_Memory& mem,Fractal& frac,vector <Point*>& hypre_points,const int& level)
  {
    //    double time0=mem.p_mess->Clock();
    const int FractalRank=mem.p_mess->FractalRank;
    const int FractalNodes=mem.p_mess->FractalNodes;
    FILE* PFH=mem.p_file->PFHypre;
    vector <int>pos(3);
    vector <int>HBox=mem.HRBoxesLev[FractalRank][level];
    vector <int>HRBox=HBox;
    vector <int>HSBox=mem.HSBoxesLev[FractalRank][level];
    fprintf(PFH," HBox %d %d %d %d %d %d \n",HBox[0],HBox[1],HBox[2],HBox[3],HBox[4],HBox[5]);
    fprintf(PFH," HSBox %d %d %d %d %d %d \n",HSBox[0],HSBox[1],HSBox[2],HSBox[3],HSBox[4],HSBox[5]);
    unsigned int minsize=mem.min_hypre_group_size;
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[level].begin();
	group_itr!=mem.all_groups[level].end();group_itr++)
      {
	Group& group=**group_itr;
	if(group.list_points.size() <= minsize || group.get_buffer_group())
	  continue;
	for(vector<Point*>::const_iterator point_itr=group.list_points.begin();point_itr !=group.list_points.end();++point_itr)
	  hypre_points.push_back(*point_itr);
      }
    //    double time1=mem.p_mess->Clock();
    int count=hypre_points.size();
    mem.p_mess->IAmAHypreNode=count > 0;
    cerr << " Hypre counts " << minsize << " " << count << " " << FractalRank << "\n";
    mem.ij_counts.clear();
    mem.ij_offsets.clear();
    mem.ij_counts.push_back(count);
    mem.ij_counts.push_back(-1);
    mem.ij_offsets.push_back(0);
    mem.ij_offsets.push_back(count);
    mem.ij_countsB=mem.ij_counts;
    mem.ij_offsetsB=mem.ij_offsets;
    //    double time2=mem.p_mess->Clock();
    mem.p_mess->HypreRank=0;
    mem.p_mess->HypreNodes=1;
    int HypreRank=0;
    int totals=count;
    if(totals == 0)
      {
	cerr << " returning " << totals << "\n";
	return false;
      }
    count=mem.ij_offsets[HypreRank];
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
    cerr << " Hypre final counts " << minsize << " " << count << " " << FractalRank << "\n";
//     cerr << " search timing " << FractalRank << " ";
//     cerr << time1-time0 << " ";
//     cerr << time2-time1 << " ";
//     cerr << time3-time2 << " ";
//     cerr << time4-time3 << " ";
//     cerr << time5-time4 << " ";
//     cerr << time6-time5 << " ";
//     cerr << time7-time6 << " ";
//     cerr << time8-time7 << "\n";
    return true;
  }
}
