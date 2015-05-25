#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int hypre_load_balance(Fractal_Memory& mem,vector <Point*>points,bool& load_balance)
  {
    FILE* PFH=mem.p_file->PFHypre;
    int HypreNodes=mem.p_mess->HypreNodes;
    int HypreRank=mem.p_mess->HypreRank;
    HYPRE_Int count_max=-1;
    double count_sum0=HypreNodes;
    double count_sum1=0.0;
    double count_sum2=0.0;
    int total=0;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	double holy_handgrenade_of_antioch=mem.ij_counts[HR];
	count_sum1+=holy_handgrenade_of_antioch;
	count_sum2+=pow(holy_handgrenade_of_antioch,2);
	count_max=max(count_max,mem.ij_counts[HR]);
	total+=mem.ij_counts[HR];
	//	fprintf(PFH," Hypre counts \t %d \t %d \n",HR,mem.ij_counts[HR]);
      }

    int average=count_sum1/count_sum0;
    int nodes_eff=count_sum1*count_sum1/count_sum2;

    bool spread_even = average >= mem.hypre_max_average_load;
    bool OOM = count_max >= mem.hypre_max_node_load;

    mem.ij_countsB=mem.ij_counts;
    fprintf(PFH," Hypre on Nodes %d %d %d \n",average,count_max,nodes_eff);
    fprintf(PFH," Hypre Load Balance %d %d \n",spread_even,OOM);
    if(!mem.hypre_load_balance)
      return 0;
    if(!OOM)
      return 0;
    load_balance=true;
    vector <HYPRE_Int> countsC=mem.ij_countsB;
    bool too_many=true;
    int trySmooth=0;
    int maxload=mem.hypre_max_node_load;
    int maxload9=(maxload*9)/10;
    while(too_many && trySmooth < 40)
      {
	too_many=false;
	for(int HR=0;HR<HypreNodes;HR++)
	  {
	    if(countsC[HR] > mem.hypre_max_node_load)
	      {
		too_many=true;
		int off=(countsC[HR]-maxload9)/5;
		//		int off=countsC[HR]/20;
		mem.ij_countsB[HR]-=2*off;
		if(HR > 0)
		  mem.ij_countsB[HR-1]+=off;
		else
		  mem.ij_countsB[HR]+=off;
		if(HR < HypreNodes-1)
		  mem.ij_countsB[HR+1]+=off;
		else
		  mem.ij_countsB[HR]+=off;
	      }
	  }
	countsC=mem.ij_countsB;
	trySmooth++;
      }
    mem.ij_offsetsB[0]=0;
    fprintf(PFH," offsets balance %d \t %d \t %d \n",0,mem.ij_offsetsB[0],mem.ij_countsB[0]);
    for(int HR=1;HR<=HypreNodes;HR++)
      {
	mem.ij_offsetsB[HR]=mem.ij_offsetsB[HR-1]+mem.ij_countsB[HR-1];
	fprintf(PFH," offsets balance %d \t %d \t %d \n",HR,mem.ij_offsetsB[HR],mem.ij_countsB[HR]);
      }
    int first_on_new_node=mem.ij_offsetsB[HypreRank];
    int last_on_new_node=mem.ij_offsetsB[HypreRank+1]-1;
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
