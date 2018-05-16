#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void move_small_boxes(Fractal_Memory& mem,vector<int>& Boxes,vector<vector<Point*>>& SPoints,vector<int>& HRout)
  {
    if(1)
      return;
    ofstream& FHT=mem.p_file->DUMPS;
    int FractalRank=mem.p_mess->FractalRank;
    int HypreRank=mem.p_mess->HypreRank;
    int HypreNodes=mem.p_mess->HypreNodes;
    const int maxLOADb=200;
    int total_boxes=accumulate(Boxes.begin(),Boxes.end(),0);
    int average_boxes=total_boxes/HypreNodes;
    int threshold=max(maxLOADb,2*average_boxes);
    int boxes_over=SPoints.size()-threshold;
    FHT << " SMALL A " << mem.level << " " << mem.steps << " " << FractalRank;
    FHT << " " << maxLOADb << " " << average_boxes << " " << threshold;
    FHT << " " << SPoints.size() << " " << boxes_over << "\n";
    if(boxes_over <= 0)
      return;
    
    // static unsigned seed1=HypreRank;
    // std::default_random_engine generator(seed1);
    // std::uniform_int_distribution<int> distribution(0,HypreNodes-1);

    multimap<int,int> boxes_here;
    for(int ni=0;ni<SPoints.size();ni++)
      boxes_here.insert(make_pair(SPoints[ni].size(),ni));
    int count=0;
    for(auto bh : boxes_here)
      {
	if(count >= boxes_over)
	  return;
	int ni=bh.second;
	cerr << " small b " << FractalRank << " " << HypreRank << " " << mem.level << " " << mem.steps << " " << boxes_over << " " << bh.first << " " << bh.second << " " << count << "\n";
	int H=HRout[ni];
	if(HRout[ni] == HypreRank)
	  HRout[ni]=rand() % HypreNodes;
	cerr << " small B " << FractalRank << " " << HypreRank << " " << mem.level << " " << mem.steps << " ";
	cerr << H << " " << HRout[ni] << "\n";
	// HRout[ni]=distribution(generator);
	count++;
      }
  }
}
