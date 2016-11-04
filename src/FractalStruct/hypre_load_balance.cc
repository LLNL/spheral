#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int hypre_load_balance(Fractal_Memory& mem,HypHest& HYP,vector <Point*>points,bool& load_balance)
  {
    FILE* PFH=mem.p_file->PFHypre;
    int HypreNodes=mem.p_mess->HypreNodes;
    int FractalRank=mem.p_mess->FractalRank;
    int HypreRank=mem.p_mess->HypreRank;
    int count_max=-1;
    double count_sum0=HypreNodes;
    double count_sum1=0.0;
    double count_sum2=0.0;
    int total=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	double holy_handgrenade_of_antioch=HYP.ij_counts[HR];
	count_sum1+=holy_handgrenade_of_antioch;
	count_sum2+=pow(holy_handgrenade_of_antioch,2);
	count_max=max(count_max,HYP.ij_counts[HR]);
	total+=HYP.ij_counts[HR];
	//	fprintf(PFH," Hypre counts \t %d \t %d \n",HR,HYP.ij_counts[HR]);
      }

    double average=count_sum1/count_sum0;
    int nodes_eff=count_sum1*count_sum1/count_sum2;
    int iaverage=average;
    int max_on_node=max((int)(average*mem.hypre_multiplier),mem.hypre_max_node_load);
//     bool spread_even = average >= max_on_node;
    bool OOM = count_max > max_on_node;
//     cerr << " loading " << FractalRank << " " << HypreRank << " " << average << " " << count_max << " " << nodes_eff << " ";
//     cerr << mem.hypre_max_node_load << " " << mem.hypre_multiplier << " " << max_on_node << " " << OOM << endl;
    HYP.ij_countsB=HYP.ij_counts;
//     fprintf(PFH," Hypre on Nodes %d %d %d %d %d %d \n",iaverage,count_max,nodes_eff,mem.hypre_max_node_load,mem.hypre_multiplier,max_on_node);
    load_balance=false;
    if(HypreNodes == 1)
      return 0;
    if(!mem.hypre_load_balance)
      return 0;
//     if(!OOM)
//       return 0;
    int loops=0;
    while(loops < 5)
      {
	for(int HR=1;HR<HypreNodes;HR++)
	  HYP.ij_offsetsB[HR]=min(HYP.ij_offsetsB[HR],HYP.ij_offsetsB[HR-1]+max_on_node);
	if(HYP.ij_offsetsB[HypreNodes]-HYP.ij_offsetsB[HypreNodes-1] <= max_on_node)
	  break;
	for(int HR=HypreNodes;HR>1;HR--)
	  HYP.ij_offsetsB[HR-1]=max(HYP.ij_offsetsB[HR-1],HYP.ij_offsetsB[HR]-max_on_node);
	if(HYP.ij_offsetsB[1]-HYP.ij_offsetsB[0] <= max_on_node)
	  break;
	loops++;
      }
    int offes=0;
    int max_off=0;
    int maxC=0;
    int maxCB=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	load_balance=load_balance || (HYP.ij_offsets[HR] != HYP.ij_offsetsB[HR]);
	int minb=min(max(HYP.ij_offsets[HR],HYP.ij_offsetsB[HR]),HYP.ij_offsets[HR+1]);
	int maxb=max(min(HYP.ij_offsets[HR+1],HYP.ij_offsetsB[HR+1]),HYP.ij_offsets[HR]);
	int offe=(minb-HYP.ij_offsets[HR])+(HYP.ij_offsets[HR+1]-maxb);
	offes+=offe;
	max_off=max(offe,max_off);
	HYP.ij_countsB[HR]=HYP.ij_offsetsB[HR+1]-HYP.ij_offsetsB[HR];
	maxC=max(maxC,HYP.ij_counts[HR]);
	maxCB=max(maxCB,HYP.ij_countsB[HR]);
 	if(HypreRank == 0)
	  fprintf(PFH," offsets balance %7d %10d %7d %10d %7d %7d %7d \n",HR,HYP.ij_offsets[HR],HYP.ij_counts[HR],HYP.ij_offsetsB[HR],HYP.ij_countsB[HR],offe,offes);
      }
    int first_on_new_node=HYP.ij_offsetsB[HypreRank];
    int last_on_new_node=HYP.ij_offsetsB[HypreRank+1]-1;
    int off_elements=0;
    for(vector<Point*>::const_iterator point_itr=points.begin();point_itr !=points.end();++point_itr)
      {
	Point* p=*point_itr;
	int label=p->get_ij_number();
	if(label < first_on_new_node || label > last_on_new_node)
	  off_elements++;
      }
    fprintf(PFH," off_elements %d %d %d %d %d %d %d \n",off_elements,max_off,offes,max_on_node,maxC,maxCB,HYP.ij_offsets[HypreNodes]);
    return off_elements;
  }
  template <class T> bool overlap_interval(T Imin,T Imax,T Jmin,T Jmax,T& LOW,T& HIGH)
  {
    HIGH=min(Imax,Jmax);
    LOW=max(Imin,Jmin);
    return HIGH >= LOW;
  }
}
