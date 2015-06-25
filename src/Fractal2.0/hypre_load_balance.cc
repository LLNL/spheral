#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int hypre_load_balance(Fractal_Memory& mem,HypHest& HYP,vector <Point*>points,bool& load_balance)
  {
    FILE* PFH=mem.p_file->PFHypre;
    int HypreNodes=mem.p_mess->HypreNodes;
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

    int average=count_sum1/count_sum0;
    int nodes_eff=count_sum1*count_sum1/count_sum2;

    bool spread_even = average >= mem.hypre_max_average_load;
    bool OOM = count_max >= mem.hypre_max_node_load;

    HYP.ij_countsB=HYP.ij_counts;
    fprintf(PFH," Hypre on Nodes %d %d %d \n",average,count_max,nodes_eff);
    fprintf(PFH," Hypre Load Balance %d %d \n",spread_even,OOM);
    if(!mem.hypre_load_balance)
      return 0;
    if(!OOM)
      return 0;
    load_balance=true;
    int trySmooth=0;
    int maxload=mem.hypre_max_node_load;
    int maxload9=(maxload*9)/10;
    int smoothMAX=(40*HypreNodes)/1024;
    smoothMAX=max(40,smoothMAX);
    bool too_many=false;
    do {
      vector <int> countsC=HYP.ij_countsB;
      too_many=false;
      for(int HR=0;HR<HypreNodes;HR++)
	{
	  if(countsC[HR] > mem.hypre_max_node_load)
	    {
	      too_many=true;
	      int off=(countsC[HR]-maxload9)/5;
	      //		int off=countsC[HR]/20;
	      HYP.ij_countsB[HR]-=2*off;
	      if(HR > 0)
		HYP.ij_countsB[HR-1]+=off;
	      else
		HYP.ij_countsB[HR]+=off;
	      if(HR < HypreNodes-1)
		HYP.ij_countsB[HR+1]+=off;
	      else
		HYP.ij_countsB[HR]+=off;
	    }
	}
      trySmooth++;
    } while(too_many && trySmooth < smoothMAX);
    HYP.ij_offsetsB[0]=0;
    fprintf(PFH," offsets balance %d \t %d \t %d \n",0,HYP.ij_offsetsB[0],HYP.ij_countsB[0]);
    for(int HR=1;HR<=HypreNodes;HR++)
      {
	HYP.ij_offsetsB[HR]=HYP.ij_offsetsB[HR-1]+HYP.ij_countsB[HR-1];
	fprintf(PFH," offsets balance %d \t %d \t %d \n",HR,HYP.ij_offsetsB[HR],HYP.ij_countsB[HR]);
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
    fprintf(PFH," off_elements %d \n",off_elements);
    return off_elements;
  }
  template <class T> bool overlap_interval(T Imin,T Imax,T Jmin,T Jmax,T& LOW,T& HIGH)
  {
    HIGH=min(Imax,Jmax);
    LOW=max(Imin,Jmin);
    return HIGH >= LOW;
  }
}
