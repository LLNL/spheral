#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_best_boxes(Fractal_Memory& mem,vector<vector<Point*> >& hypre_points,int spacing,int& VOLbest,double& FILLbest)
  {
    double time0=-mem.p_mess->Clock();
    ofstream& FHT=mem.p_file->DUMPS;
    const double maxFAKES=0.65;
    vector <double>FILLFACTOR{1.0,0.9,0.8,0.7};
    vector <int>VOLMIN {1,3,7,13,21,31,43,57,73,91,101,123,147,173,201,231};
    vector<int>MBoxes;
    vector<int>MPoints;
    vector<int>MFakes;
    for(double &FF : FILLFACTOR)
      for(int &VM : VOLMIN)
	{
	  vector<vector<int>>SBoxes;
	  vector<vector<Point*>>SPoints;
	  hypre_points_boxes(mem,hypre_points,spacing,VM,FF,SBoxes,SPoints);
	  int nPoints=0;
	  int nFakes=0;
	  for(auto &SP : SPoints)
	    {
	      for(auto &p : SP)
		{
		  nPoints++;
		  if(p == 0)
		    nFakes++;
		}
	    }
	  MBoxes.push_back(SBoxes.size());
	  MPoints.push_back(nPoints);
	  MFakes.push_back(nFakes);
	  SBoxes.clear();
	  SPoints.clear();
	}
    VOLbest=1;
    FILLbest=2.0;
    int minBOXES=999999;
    int nB=0;
    for(auto &FF : FILLFACTOR)
      for(auto &VM : VOLMIN)
	{
	  if(MBoxes[nB] < minBOXES)
	    {
	      int nFmax=(double)(MPoints[nB])*maxFAKES;
	      if(MFakes[nB] < nFmax)
		{
		  FILLbest=FF;
		  VOLbest=VM;
		  minBOXES=MBoxes[nB];
		}
	    }  
	  FHT << " BEST " << nB << " " << FF << " " << VM << " " << MBoxes[nB] << " " << MPoints[nB] << " " << MFakes[nB];
	  FHT << " " << FILLbest << " " << VOLbest << " " << minBOXES << "\n";
	  nB++;
	}
    time0+=mem.p_mess->Clock();
    FHT << " BEST TIME " << time0 << "\n";
  }
}

