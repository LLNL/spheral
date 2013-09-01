#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int hypre_load_balance(Fractal_Memory& mem,vector <Point*>points,bool& load_balance)
  {
    load_balance=false;
    ofstream& FH=mem.p_file->FileHypre;
    int FractalNodes=mem.FractalNodes;
    int FractalRank=mem.p_mess->FractalRank;
    int involved=0;
    int count_max=-1;
    double count_sum0=FractalNodes;
    double count_sum1=0.0;
    double count_sum2=0.0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	if(mem.ij_counts[FR] > 0)
	  involved++;
	double holy_handgrenade_of_antioch=mem.ij_counts[FR];
	count_sum1+=holy_handgrenade_of_antioch;
	count_sum2+=pow(holy_handgrenade_of_antioch,2);
	count_max=max(count_max,mem.ij_counts[FR]);
      }

    int average=count_sum1/count_sum0;
    int nodes_eff=count_sum1*count_sum1/count_sum2;

    bool spread_even = average >= mem.hypre_max_average_load;
    bool OOM = count_max >= mem.hypre_max_node_load;

    FH << " Hypre on Nodes " << average << " " << count_max << " " << nodes_eff << " " << involved << endl;
    if(!mem.hypre_load_balance)
      return 0;
    if(!spread_even && !OOM)
      return 0;
    FH << " Hypre Load Balance " << spread_even << OOM << endl;
    load_balance=true;
    int total=count_sum1+0.01;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	int sumb=((FR+1)*total)/FractalNodes;
	int suma=(FR*total)/FractalNodes;
	mem.ij_countsB[FR]=sumb-suma;
      }
    mem.ij_offsetsB[0]=0;
    FH << " offsets balance " << 0 << " " << mem.ij_offsets[0] << " " << mem.ij_counts[0] << endl;
    for(int FR=1;FR<=FractalNodes;FR++)
      {
	mem.ij_offsetsB[FR]=mem.ij_offsetsB[FR-1]+mem.ij_countsB[FR-1];
	FH << " offsets balance " << FR << " " << mem.ij_offsetsB[FR] << " " << mem.ij_countsB[FR] << endl;
      }
    int first_on_new_node=mem.ij_offsetsB[FractalRank];
    int last_on_new_node=mem.ij_offsetsB[FractalRank+1]-1;
    int off_elements=0;
    for(vector<Point*>::const_iterator point_itr=points.begin();point_itr !=points.end();++point_itr)
      {
	Point* p=*point_itr;
	if(p==0)
	  continue;
	int label=p->get_ij_number();
	if(label < first_on_new_node || label > last_on_new_node)
	  off_elements++;
      }
    FH << " off_elements " << off_elements << endl;
    return off_elements;
  }
  template <class T> bool overlap_interval(T Imin,T Imax,T Jmin,T Jmax,T& LOW,T& HIGH)
  {
    HIGH=min(Imax,Jmax);
    LOW=max(Imin,Jmin);
    return HIGH >= LOW;
  }
}
