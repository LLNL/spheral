#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_load_balance(Fractal_Memory& mem,vector<vector<Point*>>& SPoints)
  {
    static const double al2=1.0/log(2.0);
    ofstream& FHT=mem.p_file->DUMPS;
    int HypreRank=mem.p_mess->HypreRank;
    int HypreNodes=mem.p_mess->HypreNodes;
    vector<int>how_many(40,0);
    int maxis=-1;
    int count_points=0;
    for(auto SP : SPoints)
      {
	int is=0.01+log((double)SP.size())*al2;
	how_many[is]++;
	maxis=max(maxis,is);
	count_points+=SP.size();
      }
    vector<int>maxIS{{maxis}};
    mem.p_mess->Find_Max_INT(mem.p_mess->HypreWorld,maxIS,1);
    int MAXY=maxIS[0];
    int MAXY2=MAXY+2;
    how_many[MAXY]=SPoints.size();
    how_many[MAXY+1]=count_points;
    vector<int>boxes_on_nodes(HypreNodes*MAXY2);
    mem.p_mess->my_AllgatherI(mem.p_mess->HypreWorld,how_many,boxes_on_nodes,MAXY2);
    vector<int>Boxes(HypreNodes);
    vector<int>Points(HypreNodes);
    int ni=MAXY;
    for(int HR=0;HR<HypreNodes;HR++)
      {
	Boxes[HR]=boxes_on_nodes[ni];
	Points[HR]=boxes_on_nodes[ni+1];
	FHT << " LOAD " << HR << " ";
	for(int IS=0;IS<MAXY2;IS++)
	  FHT << boxes_on_nodes[IS+HR*MAXY2] << " ";
	FHT << "\n";
	ni+=MAXY2;
      }
    int total_boxes=accumulate(Boxes.begin(),Boxes.end(),0);
    int total_points=accumulate(Points.begin(),Points.end(),0);
    int average_boxes=total_boxes/HypreNodes;
    int average_points=total_points/HypreNodes;
    FHT << " LOADY " << total_boxes << " " << total_points << " " << average_boxes << " " << average_points << "\n";
  }
}
