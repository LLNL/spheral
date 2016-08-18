#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void box_stats(Fractal_Memory& mem,int level,int nb,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints)
  {
    ofstream& FHT=mem.p_file->DUMPS;
    int spacing=Misc::pow(2,mem.p_fractal->get_level_max()-level);
    double alog2=log(2.0);
    vector<int>counts(100,0);
    int sumP=0;
    int np=0;
    for(vector<int> &SB : SBoxes)
      {
	sumP+=SPoints[np].size();
	double VOL=1.0;
	for(int ni : {0,2,4})
	  {
	    VOL*=(double)((SB[ni+1]-SB[ni])/spacing+1);
	  }
	counts[(int)(log(VOL)/alog2)]++;
	np++;
      }
    int nimax=0;
    for(int ni=0;ni<100;ni++)
      if(counts[ni] > 0)
	nimax=ni;
    int sumFAKE=0;
    for(vector <Point*>& SP : SPoints)
      for(Point* &p : SP)
	if(p == 0)
	  sumFAKE++;
    for(int ni=0;ni<=nimax;ni++)
      {
	FHT << " BOXSTATS " << mem.p_mess->FractalRank << " " << mem.steps << " " << nb << " " << level << " " << Misc::pow(2,ni) << " " << counts[ni];
	FHT << " " << SBoxes.size() << " " << sumP << " " << sumFAKE << "\n";
      }
  }
}
