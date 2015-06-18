#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void binary_balancing(Fractal_Memory* PFM,vector <double>& numbers,double minimum,
			int Nodes,int length,vector <double>& targets,vector <int>& lowers,vector <int>& uppers)
  {
    int too_few=3;
    double VOLMAX=512.0;
//     int rank=PFM->p_mess->FractalRank;
//     double ANodes=Nodes;
    double Alength=length;
    double sum_total=std::accumulate(numbers.begin(),numbers.end(),0.0);
    vector <double>snumbers(length+1);
    minimum=1.0e-10;
    if(!PFM->periodic)
      minimum+=sum_total/Alength/(pow(VOLMAX,1.0/3.0)-1.0);
    snumbers[0]=0.0;
    for(int L=1;L<=length;L++)
      snumbers[L]=snumbers[L-1]+numbers[L-1]+minimum;
    lowers[0]=0;
    sum_total=snumbers[length];
    snumbers.resize(length);
    for(int N=1;N<Nodes;N++)
      {
	double target=sum_total*targets[N];
	lowers[N]=std::lower_bound(snumbers.begin(),snumbers.end(),target)-snumbers.begin();
	lowers[N]=max(lowers[N],lowers[N-1]+too_few);
	uppers[N-1]=lowers[N];
      }
    uppers[Nodes-1]=length-1;
    lowers[Nodes-1]=min(lowers[Nodes-1],length-too_few-1);
    uppers[Nodes-2]=lowers[Nodes-1];
    for(int N=Nodes-2;N>0;N--)
      {
	lowers[N]=min(lowers[N],lowers[N+1]-too_few);
	uppers[N-1]=lowers[N];
      }
  }
}
